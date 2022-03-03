import copy
import inspect
import sphinx.ext.autodoc as autodoc

from sphinx.ext.autodoc.mock import ismock
from sphinx.pycode import ModuleAnalyzer, PycodeError
from steps.API_2 import utils

class SubFacadeDocumenter(autodoc.ClassDocumenter):
    def __init__(self, excluded, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._excluded = excluded
        self.doc_as_attr = None

    def filter_members(self, members, want_all):
        filtered = super().filter_members(members, want_all)
        return [(name, obj, val) for name, obj, val in filtered if name not in self._excluded]

    @property
    def documenters(self):
        """Returns registered Documenter classes"""
        docs = copy.copy(self.env.app.registry.documenters)
        for name, docCls in docs.items():
            class TmpDocClass(docCls):
                def format_name(self):
                    return super().format_name().split('.')[-1]
            docs[name] = TmpDocClass
        return docs

    # Modified version of Documenter.generate
    def generate(self, more_content=None, real_modname=None,
                 check_module=False, all_members=False):
        """Generate reST for the object given by *self.name*, and possibly for
        its members.
        If *more_content* is given, include that content. If *real_modname* is
        given, use that module name to find attribute docs. If *check_module* is
        True, only generate if the object is defined in the module name it is
        imported from. If *all_members* is True, document all members.
        """
        if not self.parse_name():
            # need a module to import
            logger.warning(
                __('don\'t know which module to import for autodocumenting '
                   '%r (try placing a "module" or "currentmodule" directive '
                   'in the document, or giving an explicit module name)') %
                self.name, type='autodoc')
            return

        # now, import the module and get object to document
        if not self.import_object():
            return

        # If there is no real module defined, figure out which to use.
        # The real module is used in the module analyzer to look up the module
        # where the attribute documentation would actually be found in.
        # This is used for situations where you have a module that collects the
        # functions and classes of internal submodules.
        guess_modname = self.get_real_modname()
        self.real_modname: str = real_modname or guess_modname

        # try to also get a source code analyzer for attribute docs
        try:
            self.analyzer = ModuleAnalyzer.for_module(self.real_modname)
            # parse right now, to get PycodeErrors on parsing (results will
            # be cached anyway)
            self.analyzer.find_attr_docs()
        except PycodeError as exc:
            logger.debug('[autodoc] module analyzer failed: %s', exc)
            # no source file -- e.g. for builtin and C modules
            self.analyzer = None
            # at least add the module.__file__ as a dependency
            if hasattr(self.module, '__file__') and self.module.__file__:
                self.directive.record_dependencies.add(self.module.__file__)
        else:
            self.directive.record_dependencies.add(self.analyzer.srcname)

        if self.real_modname != guess_modname:
            # Add module to dependency list if target object is defined in other module.
            try:
                analyzer = ModuleAnalyzer.for_module(guess_modname)
                self.directive.record_dependencies.add(analyzer.srcname)
            except PycodeError:
                pass

        docstrings: List[str] = sum(self.get_doc() or [], [])
        if ismock(self.object) and not docstrings:
            logger.warning(__('A mocked object is detected: %r'),
                           self.name, type='autodoc')

        # check __module__ of object (for members not given explicitly)
        if check_module:
            if not self.check_module():
                return

        # document members, if possible
        self.document_members(all_members)


class ClassDocumenter(autodoc.ClassDocumenter):
    @staticmethod
    def getSubClasses(cls):
        subClasses = set()
        for subcls in cls.__subclasses__():
            subClasses.add(subcls)
            subClasses |= ClassDocumenter.getSubClasses(subcls)
        return subClasses

    def generate(self, *args, **kwargs):
        super().generate(*args, **kwargs)

        if issubclass(self.object, utils.Facade):
            subclasses = sorted(self.getSubClasses(self.object), key=lambda c: (inspect.getmodule(c), inspect.getsourcelines(c)[1]))
            for subcls in subclasses:
                try:
                    sectionTitle = getattr(subcls, utils.Facade._FACADE_ATTR_NAME)
                except AttributeError as e:
                    continue
                if sectionTitle is not None:
                    excluded = set()
                    for parent in subcls.__bases__:
                        excluded |= set(dir(parent))
                    remaining = set(dir(subcls)) - excluded
                    if len(remaining):
                        full_name = subcls.__module__ + '.' + subcls.__qualname__
                        self.add_line(f'**{sectionTitle}**:', self.get_sourcename())
                        self.add_line(f'   ..', self.get_sourcename())
                        documenter = SubFacadeDocumenter(excluded, self.directive, full_name, self.indent + '   ')
                        documenter.generate()


def setup(app):
    app.setup_extension('sphinx.ext.autodoc')
    app.add_autodocumenter(ClassDocumenter, override=True)

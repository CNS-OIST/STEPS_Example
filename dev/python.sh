#!/bin/bash -e


show_usage() {
    cat <<EOF
Usage: $0 [--check] [-h/--help] [path [...]]

Check Python code by calling the following utilities:

    isort: to fix sorting of module imports.
    black: to format Python code
    flake8: to perform static analysis of Python code

If invoked without option, this script will perform the following
operations on every .py of the repository:

    1. Fix import directives
    2. Format code
    3. Run static analysis

If --check option is passed, then the Python files are kept intact,
and the executable exits a non-zero value if either one of the operations
expects a change.
EOF
}

parse_options() {
    while [ "x$1" != x ] ;do
        case "$1" in
            --check)
                check=yes
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                forced_paths="${forced_paths} $1"
                ;;
        esac
        shift
    done
    if [ "x$forced_paths" != x ] ;then
        paths="$forced_paths"
    fi
}

configure() {
    eax=0
    if ! ${BLACK:-black} --version >/dev/null 2>&1 ; then
        echo "black is required. Override name with BLACK env var" >&2
        eax=$(($eax + 1))
    fi

    if ! ${ISORT:-isort} --version >/dev/null 2>&1 ; then
        echo "isort is required. Override name with ISORT env var" >&2
        eax=$(($eax + 1))
    fi

    if ! ${FLAKE8:-flake8} --version >/dev/null 2>&1 ; then
        echo "flake8 is required. Override name with FLAKE8 env var" >&2
        eax=$(($eax + 1))
    fi
    if [ $eax -ne 0 ] ;then
        echo "Abort." >&2
        exit 1
    fi
}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
paths="$SCRIPT_DIR/../python_scripts $SCRIPT_DIR/../user_manual/source"

parse_options $@
configure

eax=0

if [ "x$check" == xyes ] ;then
    if ! find $paths -name \*.py -print0 | \
        xargs -0 -n1 ${ISORT:-isort} --check-only  -w 119; then
        eax=$(($eax + 1))
    fi
    if ! ${BLACK:-black} --check -S $paths ;then
        eax=$(($eax + 1))
    fi
else
    find $paths -name \*.py -print0 | \
        xargs -0 -n1 ${ISORT:-isort} -y -w 119
    ${BLACK:-black} -S $paths
fi

if ! ${FLAKE8:-flake8} $paths ; then
    eax=$(($eax + 1))
fi

exit $eax

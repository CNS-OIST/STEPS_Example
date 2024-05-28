******************************
Customizing STEPS installation
******************************

STEPS installations can be customized by providing options to the cmake command. Some of these options were already introduced in :ref:`/getting_started.ipynb`.
The options are provided with the ``-D`` command line argument as follows::

  cmake -DOPTION_NAME=VALUE ..

CMAKE options that are specific to STEPS can be listed for your specific version with::

  cmake -LH 2> /dev/null | grep -B 1 -E '^(USE|STEPS|BUILD|ENABLE)'

Commonly useful CMAKE options
=============================

===================================== ====== ======= ===========
CMAKE option                          Type   Default Description
===================================== ====== ======= ===========
USE_MPI                               ON/OFF ON      Use MPI for parallel solvers
USE_PETSC                             ON/OFF ON      Use PETSC library for parallel E-Field solver
STEPS_INSTALL_PYTHON_DEPS             ON/OFF ON      Install the Python dependencies
STEPS_USE_DIST_MESH                   ON/OFF ON      Add solvers based on distributed mesh
STEPS_USE_HDF5_SAVING                 ON/OFF ON      Enable automatic data saving to HDF5 files
STEPS_USE_STEPSBLENDER                ON/OFF ON      Install the stepsblender python package
USE_BDSYSTEM_LAPACK                   ON/OFF OFF     Use new BDSystem/Lapack code for E-Field solver
USE_BUNDLE_EASYLOGGINGPP              ON/OFF ON      Use bundled version of easylogging
USE_BUNDLE_OMEGA_H                    ON/OFF ON      Use bundled version of Omega_h
USE_BUNDLE_RANDOM123                  ON/OFF ON      Use bundled version of random123
USE_BUNDLE_SUNDIALS                   ON/OFF ON      Use bundled version of cvode
===================================== ====== ======= ===========

CMAKE options for developpers
=============================

===================================== ====== ============ ===========
CMAKE option                          Type   Default      Description
===================================== ====== ============ ===========
BUILD_STOCHASTIC_TESTS                ON/OFF ON           Build stochastic tests
BUILD_TESTING                         ON/OFF ON           Build the testing tree
ENABLE_CODECOVERAGE                   ON/OFF OFF          Enable code coverage testing support
STEPS_ENABLE_ERROR_ON_WARNING         ON/OFF OFF          Add -Werror to STEPS compilation
STEPS_GIT_COMMIT_HOOKS                STRING ""           Comma-separated list of checks to perform when committing changes
STEPS_GIT_HOOKS                       ON/OFF ""           Enable automatic checks when committing and pushing changes
STEPS_GIT_PUSH_HOOKS                  STRING courtesy-msg Comma-separated list of checks to perform when pushing changes
STEPS_SANITIZERS                      STRING ""           Comma-separated list of runtime sanitizers to enable. Possible values: address, leak, undefined
STEPS_SANITIZERS_UNDEFINED_EXCLUSIONS STRING ""           Undefined behaviour sanitizer checks **not** to enable if STEPS_SANITIZERS contains 'undefined'
STEPS_STATIC_ANALYSIS                 ON/OFF OFF          Enable C++ static analysis during compilation
STEPS_TEST_FORMATTING                 ON/OFF OFF          Add CTest formatting test
STEPS_TEST_STATIC_ANALYSIS            ON/OFF OFF          Add CTest static analysis test
STEPS_USE_CALIPER_PROFILING           ON/OFF OFF          Use Caliper instrumentation
STEPS_USE_LIKWID_PROFILING            ON/OFF OFF          Use Likwid instrumentation
STEPS_USE_NATIVE_PROFILING            ON/OFF OFF          Use STEPS region tracker
USE_64_BITS_INDICES                   ON/OFF OFF          Use 64bits indices instead of 32
USE_CLANG_TIDY                        ON/OFF OFF          Perform C++ static analysis while compiling
===================================== ====== ============ ===========



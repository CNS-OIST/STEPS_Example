******************************
Customizing STEPS installation
******************************

STEPS installations can be customized by providing options to the cmake command. Some of these options were already introduced in :ref:`/getting_started.ipynb`.
The options are provided with the ``-D`` command line argument as follows::

  cmake -DOPTION_NAME=VALUE ..

CMAKE options that are specific to STEPS can be listed for your specific version with::

  cmake -LH 2> /dev/null | grep -B 1 -E '^(USE|STEPS)'

Commonly useful CMAKE options
=============================

===================================== ====== ===========
CMAKE option                          Type   Description
===================================== ====== ===========
USE_MPI                               ON/OFF Use MPI for parallel solvers
USE_PETSC                             ON/OFF Use PETSC library for parallel E-Field solver
STEPS_INSTALL_PYTHON_DEPS             ON/OFF Install the Python dependencies
STEPS_USE_DIST_MESH                   ON/OFF Add solvers based on distributed mesh
STEPS_USE_HDF5_SAVING                 ON/OFF Enable automatic data saving to HDF5 files
STEPS_USE_STEPSBLENDER                ON/OFF Install the stepsblender python package
USE_BDSYSTEM_LAPACK                   ON/OFF Use new BDSystem/Lapack code for E-Field solver
USE_BUNDLE_EASYLOGGINGPP              ON/OFF Use bundled version of easylogging
USE_BUNDLE_OMEGA_H                    ON/OFF Use bundled version of Omega_h
USE_BUNDLE_RANDOM123                  ON/OFF Use bundled version of random123
USE_BUNDLE_SUNDIALS                   ON/OFF Use bundled version of cvode
===================================== ====== ===========

CMAKE options for developpers
=============================

===================================== ====== ===========
CMAKE option                          Type   Description
===================================== ====== ===========
STEPS_ENABLE_ERROR_ON_WARNING         ON/OFF Add -Werror to STEPS compilation
STEPS_GIT_COMMIT_HOOKS                STRING Comma-separated list of checks to perform when committing changes
STEPS_GIT_HOOKS                       ON/OFF Enable automatic checks when committing and pushing changes
STEPS_GIT_PUSH_HOOKS                  STRING Comma-separated list of checks to perform when pushing changes
STEPS_SANITIZERS                      STRING Comma-separated list of runtime sanitizers to enable. Possible values: address, leak, undefined
STEPS_SANITIZERS_UNDEFINED_EXCLUSIONS ON/OFF Undefined behaviour sanitizer checks **not** to enable if STEPS_SANITIZERS contains 'undefined'
STEPS_STATIC_ANALYSIS                 ON/OFF Enable C++ static analysis during compilation
STEPS_TEST_FORMATTING                 ON/OFF Add CTest formatting test
STEPS_TEST_STATIC_ANALYSIS            ON/OFF Add CTest static analysis test
USE_64_BITS_INDICES                   ON/OFF Use 64bits indices instead of 32
USE_CLANG_TIDY                        ON/OFF Perform C++ static analysis while compiling
===================================== ====== ===========



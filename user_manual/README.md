# User manual

The [user manual](https://steps.sourceforge.net/manual) can be fully generated from this repository.

## Setting up

The documentation generation requires a number of python packages.
To avoid version conflicts, it is recommended to generate the documentation inside a virtual environment.

Create the virtual environment (assuming your current working directory is `user_manual`) with:
```
python3 -m venv doc_env && source doc_env/bin/activate && pip install -r requirements.txt
```

STEPS also needs to be installed in that environment. If installing from PyPI:
```
pip install STEPS
```

If installing from source, follow the normal installation procedure, but make sure that the `doc_env` environment is activated before running the `cmake ..` step.
Alternatively, one can build a source distribution locally and install it with pip (with `doc_env` activated):
```
cd /path/to/STEPS
pip install build
python -m build --sdist .
cd /path/to/STEPS_Example/user_manual
pip install /path/to/STEPS/dist/steps-5.1.0.tar.gz
```

## Generating the documentation

If needed, clean the content of the `build` directory:
```
make clean
```

Then the documentation can be generated with:
```
make html
```

## Check the documentation locally

The generated documentation can be browsed with:
```
cd build/html
python3 -m http.server
```
Then open [http://0.0.0.0:8000/](http://0.0.0.0:8000/) in a web browser.

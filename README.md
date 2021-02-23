[![Run Tests](https://github.com/HERA-Team/baseline_dependent_averaging/actions/workflows/run_tests.yaml/badge.svg)](https://github.com/HERA-Team/baseline_dependent_averaging/actions/workflows/run_tests.yaml)
[![codecov](https://codecov.io/gh/HERA-Team/baseline_dependent_averaging/branch/main/graph/badge.svg?token=JKLQhZ8mpR)](https://codecov.io/gh/HERA-Team/baseline_dependent_averaging)

# baseline_dependent_averaging
This is code for applying baseline-dependent averaging to a radio astronomy
interferometric dataset. It applies the principles and formulas presented in
[Wijnholds et
al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.2029W/abstract) to
average high-cadence data to a lower cadence while introducing a maximum amount
of decorrelation specified by the user. In brief, short baselines of an
interferometer do not decorrelate as rapidly as long baselines, and so data from
shorter baselines can be averaged together without losing as much coherent sky
information. The code and routines in this repo are designed to work on
already-recorded data (which is typically written at a common cadence for all
baselines) and averages together consecutive time samples until a specific
threshold. A forthcoming memo will describe the operation in more detail.

# Installation
The code can be installed by invoking
```sh
pip install .
```
from the top level of the repo. This will install a module called `bda` which
can be imported. The main user-facing function is `bda.apply_bda`, which is
designed to work on a [pyuvdata](https://pyuvdata.readthedocs.io/) UVData
object. It also provides a script, `apply_bda.py`, which can be called from the
command line for applying BDA to an existing dataset on disk.

## Dependencies
The following packages are required:
* astropy
* setuptools_scm
* pyuvdata

`pyuvdata` can be installed from `conda` (preferred), or from `pip`. It is
available on the `conda-forge` channel. To install:
```sh
conda install -c conda-forge pyuvdata
```

# Tests
The testing requirements can be installed by invoking
```sh
pip install .[testing]
```
from the top level of the repo. This will install the package and all
dependencies for running tests. The test suite can be run by running `pytest`
after installation.

## Dependencies
In addition to the main package dependencies above, the following packages are
required for running tests:
* pytest >= 6.0

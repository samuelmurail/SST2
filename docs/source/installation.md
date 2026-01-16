# Installation Quick Start

## Get sources from the GithubRepo

The source for `SST2` can be downloaded from the GithubRepo.

You can either clone the public repository:

```bash
$ git clone git@github.com:samuelmurail/SST2.git
```

Or download the tarball:

```bash
$ curl -OJL https://github.com/samuelmurail/SST2/tarball/master
$ tar -xvf samuelmurail-SST2-c1c2245.tar.gz
```

Once you have a copy of the source, switch to the `SST2` directory.

```bash
$ cd SST2
```

## Through Pip


Althought no pypi package is currently available, you can install `SST2` using pip:

```bash
pip install git+https://github.com/samuelmurail/SST2.git
```

## Through Conda

`SST2` can then be installed using conda:

```bash
conda env create -f environment.yml
```

or mamba:

```bash
mamba env create -f environment.yml
```

## Test installation

Use `pytest` to check that the installation was successful:

```bash
$ pip install pytest
$ pytest
=============================== test session starts ================================
platform linux -- Python 3.12.7, pytest-8.3.3, pluggy-1.5.0
rootdir: /home/murail/Documents/Code/SST2
configfile: pyproject.toml
collected 4 items                                                                  

src/SST2/tests/test_rest2.py ...                                             [ 75%]
src/SST2/tests/test_sst2.py .                                                [100%]

================================= warnings summary =================================
../../../miniforge3/envs/SST2/lib/python3.12/site-packages/pdbfixer/pdbfixer.py:58
  /home/murail/miniforge3/envs/SST2/lib/python3.12/site-packages/pdbfixer/pdbfixer.py:58: DeprecationWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html
    from pkg_resources import resource_filename

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
=========================== 4 passed, 1 warning in 37.98s =========================
```

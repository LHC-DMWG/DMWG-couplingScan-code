# DMWG-couplingScan-code
Code package implementing limit rescaling for the simplified models in arXiv:1507.00966v1. These functions allow rescaling of exclusion limits in the (mDM, mMed) plane from one coupling to another. Within certain boundaries, it also permits rescaling from one model to another. Please read [the paper](https://arxiv.org/abs/2203.12035) describing the functioning of this package in detail before using. Note that just because you can run something and get an output, it doesn't mean that output is physically meaningful!

This code expects to receive masses *in units of GeV*. If you supply masses in other units you will get nonsensical results, as the code will need to make comparisons between the mediator and dark matter masses and the masses of Standard Model particles.

## Dependencies and general installation

The package depends on the following Python packages:

```
matplotlib
pybind11
numpy
scipy
```
These dependencies should be handled automatically by the installation. 

To access the full abilities of this package (i.e. converting between models instead of only within models), you *also* need lhapdf installed. If you only want to rescale within a model and this functionality does not interest you, you can still compile without lhapdf. However, if you attempt to use any functions that rely on it, you will get errors.

Note that since pip 23.1, it is no longer supported for pip-installed packages to customise the installation with command line flags. Thus we are no longer asking for installation to specifically deactivate lhapdf. This could lead to more unexpected errors. To check if you have lhapdf available, run `lhapdf-config`. If this command exists, the installation should work. Use command `-v` to see whether lhapdf was successfully found during installation.

```
pip install -v git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
```

## Using LHAPDF on lxplus

## To install LHAPDF 

If you don't have LHAPDF available, install it. After following the section "Quick start instructions" on https://lhapdf.hepforge.org/install.html, we have to download the data files for the PDF set we want to run with. Luckily, this is quite straightforward:


```make install```
on LHAPDF instruction, do
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/folder-contains-libLHAPDF.so
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF30_nlo_as_0118.tar.gz
tar -xzvf NNPDF30_nlo_as_0118.tar.gz
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/path/to/folder-contains-NNPDF30_nlo_as_0118-folder
```

Now, do the below procedures to pip install all the relevant dependencies locally. Please do the installation in "dev mode" with lhapdf.

## Installation in development mode [instructions for Boyu and Josh]

If you want to develop

```
python -m venv thisvenv
source thisvenv/bin/activate
pip install --upgrade pip
pip install matplotlib
pip install pybind11
pip install numpy
pip install scipy
```

If you want to be able to convert vector to axial vector and vise versa, you will need lhapdf too. If you want to just test without that for now, you can install without it, but long-term it is probably helpful. If you are on lxplus, lhapdf is already available. Please edit if this is incorrect, but I believe you can use:
```
LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/
```
If it's set up correctly you shoudl also be able to run `lhapdf-config` and see some output (the help menu) indicating it correctly found that script. Unfortunately I lost lxplus access a while ago and so I can't test this ... Please feel free to make edits directly to this README once you figure out what you need to get this to work.


## Installation in development mode

If you want to modify this package locally to extend its functionality, you can install it in "dev mode". Do:

```
git clone https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
cd DMWG-couplingScan-code
pip install -e .
```

Or to install without lhapdf,
```
pip install -e . --install-option="--nolhapdf"
```

Then, see if you can successfully run the example in test/simple_test.py. 

## Usage examples

See simple working examples for different input data types in the `test` repository. These all refer to and test based on the four nominal parameter scenarios from the DMWG ( Phys.Dark Univ. 27 (2020) 100365), translating existing limits back and forth between them.

TODO paste here the text from the appendix.


## With newer pip

Required for PEP 517 and newer:

```python -m pip install -e . --config-settings="nolhapdf=True"
```
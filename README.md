# DMWG-couplingScan-code
Code package implementing limit rescaling for the simplified models in arXiv:1507.00966v1. These functions allow rescaling of exclusion limits in the (mDM, mMed) plane from one coupling to another. Within certain boundaries, it also permits rescaling from one model to another. Please read [the paper](insert link when submitted) describing the functioning of this package in detail before using. Note that just because you can run something and get an output, it doesn't mean that output is physically meaningful!

This code expects to receive masses *in units of GeV*. If you supply masses in other units you will get nonsensical results, as the code will need to make comparisons between the mediator and dark matter masses and the masses of Standard Model particles.

## Dependencies and general installation

The package depends on the following Python packages:

```
pybind11
numpy
scipy
```
These dependencies should be handled automatically by the installation. If this doesn't work (e.g. due to your cluster configuration), install them yourself before proceeding.

To access the full abilities of this package (i.e. converting between models instead of only within models), you *also* need lhapdf installed. If you only want to rescale within a model and this functionality does not interest you, you can still compile without lhapdf. However, if you attempt to use any functions that rely on it, you will get errors.

Note that since pip 23.1, it is no longer supported for pip-installed packages to customise the installation with command line flags. Thus we are no longer asking for installation to specifically deactivate lhapdf. This could lead to more unexpected errors. To check if you have lhapdf available, run `lhapdf-config`. If this command exists, the installation should work. Use command `-v` to see whether lhapdf was successfully found during installation.

```
python -m pip install -v git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
```

****


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

## Installation in development mode (for contributors to the package)

If you want to develop the code, follow this:

```
python -m venv thisvenv
source thisvenv/bin/activate
pip install --upgrade pip
pip install matplotlib
pip install pybind11
pip install numpy
pip install scipy
```

And you need lhapdf too. If you are on lxplus, lhapdf is already available. Please edit if this is incorrect, but I believe you can use:
```
LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/
```
If it's set up correctly you shoudl also be able to run `lhapdf-config` and see some output (the help menu) indicating it correctly found that script. 

If you're on your laptop and you want lhapdf, you'll have to install it yourself following the instructions above. You can use the `lhapdf-config` test again.

Now, to install in "dev mode", that is so that you get a local copy of the code and can develop it, do this:

```
git clone https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
cd DMWG-couplingScan-code
pip install -e .
```

You should be able to run all scripts in the `test` directory.

To test/use as a general user would do, you can use the pip install instructions in the previous section. 

## Usage examples

See simple working examples for different input data types in the `test` repository. These all refer to and test based on the four nominal parameter scenarios from the DMWG ( Phys.Dark Univ. 27 (2020) 100365), translating existing limits back and forth between them.

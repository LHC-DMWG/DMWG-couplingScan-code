# DMWG-couplingScan-code
Code package implementing limit rescaling for the simplified models in arXiv:1507.00966v1. These functions allow rescaling of exclusion limits in the (mDM, mMed) plane from one coupling to another. Within certain boundaries, it also permits rescaling from one model to another. Please read [the paper](https://arxiv.org/abs/2203.12035) describing the functioning of this package in detail before using. Note that just because you can run something and get an output, it doesn't mean that output is physically meaningful!

This code expects to receive masses *in units of GeV*. If you supply masses in other units you will get nonsensical results, as the code will need to make comparisons between the mediator and dark matter masses and the masses of Standard Model particles.

## Dependencies and general installation

The package depends on:

```
matplotlib
pybind11
numpy
scipy
```
You'll need to have all of these installed. 
To access the full abilities of this package (i.e. converting between models instead of only within models), you also need lhapdf. If you only want to rescale within a model and this functionality does not interest you, you can compile with the `no-lhapdf` flag:

```
pip install git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
```
or
```
pip install git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git --install-option="--no-lhapdf"
```

****

## Using LHAPDF on lxplus

If you are an lxplus user, lhapdf is already available. I think you need to do this:
```
LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/
```
If it's set up correctly you should be able to run `lhapdf-config` and see some output (the help menu) indicating it correctly found that script. Unfortunately I do not have lxplus access any longer, so I can't keep these instructions up to date. Please feel free to make edits directly to this README if you need alternate instructions to get this to work.

## Install LHAPDF manually

If you are not running on lxplus and you want full functionality, you will need to install LHAPDF on your local machine. After following the section "Quick start instructions" on https://lhapdf.hepforge.org/install.html, download the data files and link the dynamical path and data path to the right place. i.e. after
```make install```
on LHAPDF instruction, do
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/folder-contains-libLHAPDF.so
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF30_nlo_as_0118.tar.gz
tar -xzvf NNPDF30_nlo_as_0118.tar.gz
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/path/to/folder-contains-NNPDF30_nlo_as_0118-folder
```

You can use the `lhapdf-config` to test your setup.

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

See simple working examples here: https://github.com/LHC-DMWG/DMWG-couplingScan-code/blob/master/test/simple_test.py

TODO extend this section

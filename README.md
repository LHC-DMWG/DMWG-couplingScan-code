# DMWG-couplingScan-code
Code package implementing limit rescaling for the simplified models in arXiv:1507.00966v1. These functions allow rescaling of exclusion limits in the (mDM, mMed) plane from one coupling to another. Within certain boundaries, it also permits rescaling from one model to another. Please read [the paper](https://arxiv.org/abs/2203.12035) describing the functioning of this package in detail before using. Note that just because you can run something and get an output, it doesn't mean that output is physically meaningful!

This code expects to receive masses *in units of GeV*. If you supply masses in other units you will get nonsensical results, as the code will need to make comparisons between the mediator and dark matter masses and the masses of Standard Model particles.

## Dependencies and general installation

The package depends on:

```
pybind11
numpy
scipy
```

To access the full abilities of this package (i.e. converting between models instead of only within models), you need lhapdf installed. If you only want to rescale within a model and this functionality does not interest you, you can compile with the `no-lhapdf` flag:

```
pip install git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
```
or
```
pip install git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git --install-option="--no-lhapdf"
```

****


## Install LHAPDF manually [another way of installation]

After following the section "Quick start instructions" on https://lhapdf.hepforge.org/install.html, we have to download the data files and link the dynamical path and data path to the right place. i.e. after
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

Hi guys! You can do the above if you want but you will probably prefer to have access to the code locally. I would suggest you do this:

```
python -m venv thisvenv
source thisvenv/bin/activate
pip install --upgrade pip
pip install matplotlib
pip install pybind11
pip install numpy
pip install scipy
```

If you want to be able to convert vector to axial vector and vise versa, you will need lhapdf too. If you want to just test without that for now, you can install without it, but long-term it is probably helpful. If you are on lxplus, lhapdf is already available. If i remember correctly, you need to do this:
```
LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/
```
If it's set up correctly you shoudl also be able to run `lhapdf-config` and see some output (the help menu) indicating it correctly found that script. Unfortunately I lost lxplus access a while ago and so I can't test this ... Please feel free to make edits directly to this README once you figure out what you need to get this to work.

If you're on your laptop and you want lhapdf, you'll have to install it yourself. You can use the `lhapdf-config` test again.

Now, to install in "dev mode", that is so that you get a local copy of the code and can develop it, do this:

```
git clone https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
cd DMWG-couplingScan-code
pip install -e .
```

Or to install without lhapdf,
```
pip install -e . --install-option="--nolhapdf"
```

Then, see if you can successfully run the example in test/simple_test.py. These numbers aren't validated yet - I am hoping you can both help with that! - but all the tests should run successfully and give real number outputs if things are working well. If you don't have lhapdf installed, one or two of these tests will fail, but the error message will be transparent.

To test/use as a general user would do, you can use the pip install instructions in the previous section. 

## Usage examples

See simple working examples here: https://github.com/LHC-DMWG/DMWG-couplingScan-code/blob/master/test/simple_test.py

TODO extend this section

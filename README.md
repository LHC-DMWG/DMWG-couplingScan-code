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

To access the full abilities of this package (i.e. converting between models instead of only within models), you *also* need lhapdf installed. If you only want to rescale within a model and this functionality does not interest you, you can still compile without lhapdf. However, if you attempt to use any functions that rely on it, the program will throw an error and exit. Note that since pip 23.1, it is no longer supported for pip-installed packages to customise the installation with command line flags. Thus we are no longer asking for the user to specifically deactivate lhapdf when installing. 

To check if you have lhapdf available, run `lhapdf-config`. If this command exists, the installation should work.

```
python -m pip install git+https://github.com/LHC-DMWG/DMWG-couplingScan-code.git
```

## Installing/accessing LHAPDF 

If you don't have LHAPDF available, install it. Follow the "Quick start instructions" on https://lhapdf.hepforge.org/install.html. If you are on lxplus, lhapdf is already available. Please edit if this instruction is incorrect, but I believe you can use:
```
LHAPDF_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/
''' 

You may also have to download the data files for the PDF set you want to run with (by default, NNPDF30_nlo_as_0118). Luckily, this is quite straightforward using the tools provided:

```lhapdf get NNPDF30_nlo_as_0118```
adding where to download (`--pdfdir`) as required. You can alternatively download directly as a tarball:

```
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF30_nlo_as_0118.tar.gz
tar -xzvf NNPDF30_nlo_as_0118.tar.gz
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/path/to/folder-containsing-NNPDF30_nlo_as_0118
```

## When running

Integration warnings are fairly common when running the mono-X scans, with two or three often arising per grid studied. These do not so far appear to cause any real numerical instabilities in the results or have visible consequences, so 

## Usage examples

See simple working examples for different input data types in the `test` repository. These all refer to and test based on the four nominal parameter scenarios from the DMWG (Phys.Dark Univ. 27 (2020) 100365), translating existing limits back and forth between them.

TODO paste here the text from the appendix.

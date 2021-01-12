from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
import os

# Run this script with:
# python setup.monox_functions.py build_ext --inplace

lhapdf_path = os.environ.get('LHAPDF_LIBRARY_PATH')
print("path is",lhapdf_path)

setup(
    ext_modules=cythonize(Extension("monox_functions",
                                    sources=["monox_functions.pyx"],
                                    include_dirs=['./','/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0-cebb0/x86_64-centos7/include/c++/8.3.0','/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0-cebb0/x86_64-centos7/include/c++/8.3.0/x86_64-pc-linux-gnu'],
                                    libraries=["LHAPDF"],
                                    library_dirs=["/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/lib"]
                                    ),
                          annotate=True)
)


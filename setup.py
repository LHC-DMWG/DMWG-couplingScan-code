from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Detailed examples: 
# https://github.com/pybind/pybind11_benchmark/blob/master/setup.py
# https://github.com/wichert/pybind11-example/blob/master/setup.py
def get_lhapdf_includes() :
   import subprocess
   linker_flags = subprocess.check_output(["lhapdf-config", "--libs"])
   print(linker_flags)

#get_lhapdf_includes()

__version__ = "0.0.1"
ext_modules = [
    Pybind11Extension(
        "lhapdfwrap",
        sources = ["src/lhapdf_integrands.cpp"],
        library_dirs = ['/usr/local/Cellar/lhapdf/6.2.1/lib'],
        libraries = ['LHAPDF'],
        include_dirs = ['/usr/local/Cellar/lhapdf/6.2.1/include'],
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name = 'package',
    version=__version__,
    url = 'https://github.com/LHC-DMWG/DMWG-couplingScan-code/',
    author = 'LHC DMWG',
    author_email = 'lhc-dmwg-contributors-couplingscan@cern.ch',
    description = '',
    packages = find_packages(),    
    install_requires = requirements,
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext}    
)
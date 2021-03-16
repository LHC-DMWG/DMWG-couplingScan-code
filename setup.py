from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Detailed examples: 
# https://github.com/pybind/pybind11_benchmark/blob/master/setup.py
# https://github.com/wichert/pybind11-example/blob/master/setup.py
def get_lhapdf_includes() :
    libdir = subprocess.check_output(["lhapdf-config", "--libdir"])
    incdir = subprocess.check_output(["lhapdf-config", "--libdir"])
    if not libdir or not incdir :
        print("You must have lhapdf installed to use this tool!")
        exit(1)
    return libdir,incdir


lhapdf_dirs = get_lhapdf_includes()

__version__ = "0.0.1"
ext_modules = [
    Pybind11Extension(
        "lhapdfwrap",
        sources = ["src/lhapdf_integrands.cpp"],
        library_dirs = [lhapdf_dirs[0]],
        libraries = ['LHAPDF'],
        include_dirs = [lhapdf_dirs[1]],
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
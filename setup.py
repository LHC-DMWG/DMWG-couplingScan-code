from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

__version__ = "0.0.1"
ext_modules = [
    Pybind11Extension("lhapdfwrap",
        ["src/main.cpp"],
        # Example: passing in the version to the compiled code
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
    ext_modules=ext_modules
)

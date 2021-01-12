from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name = 'package',
    version = '0.0.1',
    url = 'https://github.com/LHC-DMWG/DMWG-couplingScan-code/',
    author = 'LHC DMWG',
    author_email = 'lhc-dmwg-contributors-couplingscan@cern.ch',
    description = '',
    packages = find_packages(),    
    install_requires = requirements,
)

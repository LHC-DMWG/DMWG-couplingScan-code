from distutils.core import setup
from distutils.extension import Extension

import subprocess
def find_lhapdf() :
    command = "lhapdf-config --"

setup(name="Test",
      ext_modules=[
          Extension("monox_functions", 
                    sources = ["monox_functions.cpp"],
                    libraries = )

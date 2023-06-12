from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import subprocess
import sys

# Detailed examples: 
# https://github.com/pybind/pybind11_benchmark/blob/master/setup.py
# https://github.com/wichert/pybind11-example/blob/master/setup.py
def found_lhapdf() :
    try :
        config = subprocess.check_output(["lhapdf-config"]).decode('ascii').strip()
        print("Got config:")
        print(config)
        return True
    except :
        config = ""
    return config != ""

def get_lhapdf_includes() :
        
    libdir = subprocess.check_output(["lhapdf-config", "--libdir"]).decode('ascii').strip()
    incdir = subprocess.check_output(["lhapdf-config", "--incdir"]).decode('ascii').strip()

    return libdir,incdir

def get_extmodules() :

    # Check if we are actually compiling the C++ code.
    lhapdf = found_lhapdf()
    if not lhapdf :
        ext_modules = []
        print("""\nRequested installation without LHAPDF dependency. 
        Note that this will limit the functionality!\n""")
    else :
        # Here I need to access the flags.
        lhapdf_dirs = get_lhapdf_includes()
        ext_modules = [
            Pybind11Extension(
                "lhapdfwrap",
                sources = ["src/lhapdf_integrands.cpp"],
                library_dirs = [lhapdf_dirs[0]],
                libraries = ['LHAPDF'],
                include_dirs = [lhapdf_dirs[1]],
                ),
        ]
    return ext_modules

setup(
    ext_modules=get_extmodules()
)

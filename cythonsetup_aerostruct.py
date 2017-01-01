from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
	ext_modules = cythonize("aerostruct.pyx"),
	include_dirs=[numpy.get_include()]
)

# Run this from the Anaconda prompt to compile Cython code...
# $ python cythonsetup_aerostruct.py build_ext --inplace

from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
	ext_modules = cythonize("cytn/aerostruct_cython.pyx"),
	include_dirs=[numpy.get_include()]
)

# Run this from the Anaconda prompt to compile Cython code...
# $ python setup_aerostruct_cython.py build_ext --inplace

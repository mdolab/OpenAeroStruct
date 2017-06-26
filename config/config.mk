# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

# Fortran compiler and flags
FF90        = gfortran
FF90_FLAGS  = -fdefault-real-8 -O2 -fPIC

# C compiler and flags
CC       = gcc
CC_FLAGS   = -O2 -fPIC

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS = -fPIC -shared

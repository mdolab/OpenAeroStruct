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
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
ifeq ($(python_version_major), 2)
    PYTHON-CONFIG = python-config
endif
ifeq ($(python_version_major), 3)
    PYTHON-CONFIG = python3.6-config
endif
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS = -fPIC -shared

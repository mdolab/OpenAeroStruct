# ----------------------------------------------------------------------
# Config file for Gfortran  with OpenMPI
# ----------------------------------------------------------------------

# ------- Define a possible parallel make ------------------------------
PMAKE = make -j 4

# ------- Define the MPI Compilers--------------------------------------
FF90 = gfortran
CC   = cc

# ------- Define Compiler Flags ----------------------------------------
FF90_GEN_FLAGS = -fPIC -g -fbounds-check
CC_GEN_FLAGS   = -fPIC

FF90_OPT_FLAGS   = -fPIC -fdefault-real-8 -O2 -fdefault-double-8
CC_OPT_FLAGS     = -O2

# ------- Define Archiver  and Flags -----------------------------------
AR       = ar
AR_FLAGS = -rvs

# ------- Define Linker Flags ------------------------------------------
LINKER_FLAGS =

# Use the following line for Mac install
# LINKER_FLAGS = -undefined dynamic_lookup

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python-config
F2PY = f2py

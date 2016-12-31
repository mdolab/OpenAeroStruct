#/bin/sh

SRC_FILES = src/aeromtx.f90 \
					src/structmtx.f90

default:
	python -m numpy.f2py -c ${SRC_FILES} -m lib

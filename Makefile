#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author: John Jasa                  							              *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 12-08-2016                                      *
#      *                                                                *
#      ******************************************************************

SUBDIR_SRC    =	src/adjoint \
		src/OAS \

default:

# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_INTEL_OPENMPI.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. Typically the CGNS directory "; \
	echo "will have to be modified. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make oas;\
	fi;

clean:
	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f mod/* obj/*
	(rm *.so) || exit 1;

oas:
	mkdir -p obj
	mkdir -p mod
	ln -sf config/config.mk config.mk
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)
	(cd src/python/f2py && make)
	mkdir -p python
	(cd python && cp *.so ../)

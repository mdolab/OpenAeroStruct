# Makefile for pySpline

SUBDIR_SRC    = src/adjoint \
							  src/OAS

default:

# Check if the config.mk file is in the config dir.
	@if [ ! -f "config/config.mk" ]; then \
	echo "Before compiling, copy an existing config file from the "; \
	echo "config/defaults/ directory to the config/ directory and  "; \
	echo "rename to config.mk. For example:"; \
	echo " ";\
	echo "  cp config/defaults/config.LINUX_INTEL.mk config/config.mk"; \
	echo " ";\
	echo "The modify this config file as required. With the config file specified, rerun "; \
	echo "'make' and the build will start"; \
	else make module;\
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
	rm -f lib/lib* mod/*.mod obj/*

module:
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

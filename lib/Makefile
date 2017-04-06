#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: John Jasa          							                      *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 07-27-2016                                      *
#      *                                                                *
#      ******************************************************************
#
HOME_DIR = ..
COMMON_FILE = $(HOME_DIR)/common.mk
include ${COMMON_FILE}

default: all

all: 	liboas

clean:
	@echo "        Making clean in lib ... "
	rm -f lib/*

liboas:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) liboas.a $(OBJDIR)/*.o
	@echo
	@echo "        library liboas.a created."
	@echo

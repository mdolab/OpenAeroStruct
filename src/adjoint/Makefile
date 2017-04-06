#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Author: John Jasa                 					  		              *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 12-05-2016                                      *
#      *                                                                *
#      ******************************************************************

HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = adBuffer.o \
		 adStack.o \
		 outputReverse/oas_main_b.o \
		 outputForward/oas_main_d.o

default: all

all:
	@echo "            Making in OpenAeroStruct Derivatives... "
	# Only possible if you have Tapenade installed
	# $(PMAKE)  derivatives
	$(PMAKE)  src1
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in OpenAeroStruct Derivatives... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

derivatives:
	@echo "     Running Tapenade in OpenAeroStruct Derivatives... "
	# Make the forward mode derivative files
	make -f Makefile_tapenade

	# Make the reverse mode derivative files
	make -f Makefile_tapenade ad_reverse

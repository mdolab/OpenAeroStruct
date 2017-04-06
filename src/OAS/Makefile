#      ******************************************************************
#      *                                                                *
#      * File:          Makefile                                        *
#      * Authors: John Jasa                      							          *
#      * Based on Gaetan Kenway's Makefiles                             *
#      * Starting date: 07-27-2016                                      *
#      * Last modified: 10-31-2016                                      *
#      *                                                                *
#      ******************************************************************

HOME_DIR = ../..
COMMON_FILE = $(HOME_DIR)/common.mk
RULES_FILE  = $(HOME_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = OAS_main.o \

FF90_OBJECTS_2 = OAS_API.o

default: all

all:
	$(PMAKE)  src1
	$(PMAKE)  src2
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in OpenAeroStruct... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

src1:	$(FF90_OBJECTS_1)

src2:	$(FF90_OBJECTS_2)

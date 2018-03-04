MAIN_DIR = ../..
COMMON_FILE = $(MAIN_DIR)/common.mk
RULES_FILE  = $(MAIN_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = adBuffer.o \
		 adStack.o \
		 outputReverse/oas_main_b.o \
		 outputForward/oas_main_d.o

default:
	@echo "            Making in OpenAeroStruct Derivatives... "
	make  $(FF90_OBJECTS_1)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in OpenAeroStruct Derivatives... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

derivatives:
	@echo "     Running Tapenade in OpenAeroStruct Derivatives... "
	# Make the forward mode derivative files
	make -f Makefile_tapenade

	# Make the reverse mode derivative files
	make -f Makefile_tapenade ad_reverse

MAIN_DIR = ../..
COMMON_FILE = $(MAIN_DIR)/common.mk
RULES_FILE  = $(MAIN_DIR)/rulesSources.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

FF90_OBJECTS_1 = OAS_main.o \

FF90_OBJECTS_2 = OAS_API.o

default:
	make  $(FF90_OBJECTS_1)
	make  $(FF90_OBJECTS_2)
	cp -f *.mod $(MODDIR)

clean:
	@echo "        Making clean in OpenAeroStruct... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

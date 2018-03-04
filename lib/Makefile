MAIN_DIR = ..
COMMON_FILE = $(MAIN_DIR)/common.mk
include ${COMMON_FILE}

default: makelib

clean:
	@echo "        Making clean in lib ... "
	rm -f lib/*

#      ******************************************************************
#      *                                                                *
#      * Rule to make the library.                                      *
#      *                                                                *
#      ******************************************************************

makelib:
	@echo "        creating library ..."
	@echo
	$(AR) $(AR_FLAGS) liboas.a $(OBJDIR)/*.o
	@echo
	@echo "        library liboas.a created."
	@echo

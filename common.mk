MODDIR = $(HOME_DIR)/mod
OBJDIR = $(HOME_DIR)/obj
LIBDIR = $(HOME_DIR)/lib

#      ******************************************************************
#      *                                                                *
#      * Include the file describing the compiler settings.             *
#      *                                                                *
#      ******************************************************************

COMPILERS = $(HOME_DIR)/config.mk
ifneq ($(MAKECMDGOALS),clean)
include ${COMPILERS}
endif

#      ******************************************************************
#      *                                                                *
#      * Redefine .SUFFIXES to be sure all the desired ones are         *
#      * included.                                                      *
#      *                                                                *
#      ******************************************************************

.SUFFIXES: .o .f .F .f90 .F90

#      ******************************************************************
#      *                                                                *
#      * Arguments of make clean.                                       *
#      *                                                                *
#      ******************************************************************

MAKE_CLEAN_ARGUMENTS = *~ *.o *.mod *.il *.stb c_* *.a *.so

#      ******************************************************************
#      *                                                                *
#      * Compiler flags to compile the sources.                         *
#      * The macro's ADDITIONAL_FF90_FLAGS and ADDITIONAL_CC_ALL_FLAGS  *
#      * make it possible that every subdirectory adds its own specific *
#      * compiler flags, if necessary.                                  *
#      *                                                                *
#      ******************************************************************

FF90_ALL_FLAGS   = -I$(MODDIR) $(FF90_GEN_FLAGS) $(FF90_OPT_FLAGS)

CC_ALL_FLAGS = -I$(MODDIR) $(CC_GEN_FLAGS) $(CC_OPT_FLAGS)

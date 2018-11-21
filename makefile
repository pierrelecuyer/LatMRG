# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++14 -Wall -O2
DEBUG_FLAGS = -std=c++14 -g -Wall -O2

# The header files of the LatMRG library and the LatticeTester library
INCLUDES = -I./include -I./latticetester/include

# This is included for backwards compatibility
DEF_LLDD = -DNTL_TYPES_CODE=1
DEF_ZZDD = -DNTL_TYPES_CODE=2
DEF_ZZRR = -DNTL_TYPES_CODE=3
NUM_TYPES = $(DEF_ZZDD)

# Library path. This assumes NTL is in /usr/local/lib (its default path).
STAT_LIBS_PATH = -Wl,-Bstatic -L$(LIB_DIR)
STAT_LIBS = -llatmrg -llatticetester 
DYN_LIBS_PATH = -Wl,-Bdynamic -L/usr/local/lib 
DYN_LIBS = -lntl -lgmp #-ltestu01 -lmylib

# A few directories we need to be aware of
SRC_DIR = ./src
OBJ_DIR = ./obj
LIB_DIR = ./lib
BIN_DIR = ./bin
INC_DIR = ./include
PRO_DIR = ./progs

EX_DIR = ./examples
EX_BUILD = ./bin/examples

# Other source files locations
MRG_LOW = $(SRC_DIR)/latmrg
MRG_HIGH  = $(SRC_DIR)/latmrg-high
MRG_TYPES = $(SRC_DIR)/latmrg/mrgtypes

# The source files are in SRC_DIR. This grabs subdirectories
SRCS = $(wildcard $(MRG_LOW)/*.cc) $(wildcard $(MRG_HIGH)/*.cc) \
       $(wildcard $(MRG_TYPES)/*.cc)
PROGS_CC = $(wildcard $(PRO_DIR)/*.cc)
EX_CC = $(wildcard $(EX_DIR)/*.cc)

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
PROGS_O = $(PROGS_CC:$(PRO_DIR)/%.cc=$(PRO_DIR)/%.o)
EX_O = $(EX_CC:%.cc=%.o)

# A separator to segment the information printed on screen

SEP = @echo ================================================================================

default: lib progs

all: clean_all lib bin examples doc

bin: lib progs

#==============================================================================
# Building the API

lib: latticetester $(LIB_DIR)/ lib_objects
	cp latticetester/build/src/liblatticetester.a $(LIB_DIR)
	rm -f $(LIB_DIR)/liblatmrg.a
	ar rcs $(LIB_DIR)/liblatmrg.a $(OBJS)
	@echo 'LatMRG library archive created in ./lib/liblatmrg.a'
	@echo

$(LIB_DIR)/:
	$(SEP)
	@echo 'Building the LatMRG library'
	@echo
	mkdir -p $(LIB_DIR)

$(OBJ_DIR)/:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/latmrg
	mkdir -p $(OBJ_DIR)/latmrg/mrgtypes
	mkdir -p $(OBJ_DIR)/latmrg-high

lib_objects: $(OBJ_DIR)/ $(OBJS)

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) $(NUM_TYPES) -c $< -o $@

#==============================================================================
# Building the programs of ./progs

progs: $(BIN_DIR)/ progs_objects
	@echo 'LatMRG programs compiled in ./bin'
	@echo

$(BIN_DIR)/:
	$(SEP)
	@echo 'Building the LatMRG executables'
	@echo
	mkdir -p $(BIN_DIR)

# Builds objects and binaries for every of the programs of LatMRG
progs_objects: $(PROGS_O)

$(PRO_DIR)/%.o:$(PRO_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(NUM_TYPES) -c $< -o $@
	$(CC) $@ $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) $(DYN_LIBS) \
	  -o $(BIN_DIR)/$(@:progs/%.o=%) 

#==============================================================================
# Building the documentation

doc:
	$(SEP)
	@echo 'Building the documentation'
	@echo
	doxygen doc/doc-gen
	@echo 'Documentation now in ./doc'
	@echo

#==============================================================================
# Building the examples

examples:$(EX_BUILD)/ build_ex
	@echo 'LatMRG examples compiled in ./bin/examples'
	@echo

$(EX_BUILD)/:
	@echo
	$(SEP)
	@echo 'Building the examples'
	@echo
	mkdir -p $(EX_BUILD)

build_ex:$(EX_O)

$(EX_DIR)/%.o:$(EX_DIR)/%.cc
	$(CC) $< $(INCLUDES) -I. $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) \
	  $(DYN_LIBS) -o $(EX_BUILD)/$(<:examples/%.cc=%) 

#==============================================================================
# Installation/removal of LatMRG 

#==============================================================================
# LatticeTester related considerations

config_latticetester:
	@read -p 'Enter NTL prefix (default: /usr/local): ' NTL_PREFIX; \
	  if [ -z $$NTL_PREFIX ]; then\
	    NTL_PREFIX='--ntl /usr/local';\
	  else\
	    NTL_PREFIX='--ntl '$$NTL_PREFIX;\
	  fi;\
	  read -p 'Enter boost prefix (default: empty string): ' BOOST_PREFIX;\
	  if [ -z $$BOOST_PREFIX ] ; then\
	    BOOST_PREFIX='';\
	  else\
	    BOOST_PREFIX='--boost '$$BOOST_PREFIX;\
	  fi;\
	  cd latticetester;\
	  ./waf configure $$NTL_PREFIX $$BOOST_PREFIX

latticetester:
	$(SEP)
	@echo 'Building LatticeTester with waf'
	@echo
	cd latticetester; ./waf build
	@echo 'LatticeTester build finished'
	@echo

clean_lattice:
	$(SEP)
	@echo 'Cleaning LatticeTester'
	@echo
	cd latticetester; ./waf clean
	@echo

#==============================================================================
# Cleaning stuff

clean: separator clean_lib clean_progs clean_doc
	@echo

clean_all: clean_lattice clean

clean_lib:
	rm -rf $(LIB_DIR)
	rm -rf $(OBJ_DIR)

clean_progs:
	rm -rf $(BIN_DIR)
	rm -f $(PRO_DIR)/*.o

clean_doc:
	rm -rf ./doc/html
	rm -rf ./doc/latex

clean_examples:
	rm -rf $(BIN_DIR)/examples

#==============================================================================
# Graphical stuff

separator:
	$(SEP)
	@echo

#==============================================================================
# PHONY targets

.PHONY: latticetester $(LIB_DIR)/ $(OBJ_DIR)/ $(BIN_DIR)/ $(OBJ_DIR)/ doc clean\
  clean_all examples $(EX_BUILD)/ $(EX_CC) separator config_latticetester

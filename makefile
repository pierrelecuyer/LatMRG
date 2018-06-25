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

default: lib progs

all: clean_all lib bin examples doc

bin: lib progs

#==============================================================================
# Building the API

lib: latticetester $(LIB_DIR)/ lib_objects
	cp latticetester/lib/liblatticetester.a $(LIB_DIR)
	rm -f $(LIB_DIR)/liblatmrg.a
	ar rcs $(LIB_DIR)/liblatmrg.a $(OBJS)

$(LIB_DIR)/:
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

$(BIN_DIR)/:
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
	doxygen doc/doc-gen

#==============================================================================
# Building the examples

examples:$(EX_BUILD)/ build_ex

$(EX_BUILD)/:
	mkdir -p $(EX_BUILD)

build_ex:$(EX_O)

$(EX_DIR)/%.o:$(EX_DIR)/%.cc
	$(CC) $< $(INCLUDES) $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) \
	  $(DYN_LIBS) -o $(EX_BUILD)/$(<:examples/%.cc=%) 

#==============================================================================
# Installation/removal of LatMRG 

#==============================================================================
# LatticeTester related considerations

latticetester:
	cd latticetester; make

clean_lattice:
	cd latticetester; make clean

#==============================================================================
# Cleaning stuff

clean: clean_lib clean_progs clean_doc

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
# PHONY targets

.PHONY: latticetester $(LIB_DIR)/ $(OBJ_DIR)/ $(BIN_DIR)/ $(OBJ_DIR)/ doc clean\
  clean_all examples $(EX_BUILD)/ $(EX_CC)

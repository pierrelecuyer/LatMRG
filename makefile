# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++14 -Wall -O2
DEBUG_FLAGS = -std=c++14 -g -Wall -O2

# define any directories containing header files other than /usr/include
INCLUDES = -I./include -I./latticetester/include

# This is included for backwards compatibility
DEF_LLDD = -DNTL_TYPES_CODE=1
DEF_ZZDD = -DNTL_TYPES_CODE=2
DEF_ZZRR = -DNTL_TYPES_CODE=3
NUM_TYPES = $(DEF_ZZDD)

# define library paths in addition to /usr/lib
STAT_LIBS_PATH = -Wl,-Bstatic -L$(LIB_DIR)
STAT_LIBS = -llatmrg -llatticetester 
DYN_LIBS_PATH = -Wl,-Bdynamic -L/usr/local/lib 
DYN_LIBS = -lntl -lgmp #-ltestu01 -lmylib

# A few directories we need to be aware of
SRC_DIR = ./src
OBJ_DIR = ./obj
LIB_DIR = ./lib
BIN_DIR = ./bin
PRO_DIR = ./progs


# define the C source files
SRCS = $(wildcard $(SRC_DIR)/*.cc)
PROGS_CC = $(wildcard $(PRO_DIR)/*.cc)

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

all: mkdir lib bin

lib: latticetester lib_objects
	rm -f $(LIB_DIR)/liblatmrg.a
	ar rcs $(LIB_DIR)/liblatmrg.a $(OBJS)

lib_objects: $(OBJS)

bin: lib progs_objects

progs_objects: $(PROGS_O)

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(NUM_TYPES) -c $< -o $@

$(PRO_DIR)/%.o:$(PRO_DIR)/%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(NUM_TYPES) -c $< -o $@
	$(CC) $@ $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) $(DYN_LIBS) -o $(BIN_DIR)/$(@:progs/%.o=%) 

#==============================================================================

latticetester:
	cd latticetester; make
	cp latticetester/lib/liblatticetester.a $(LIB_DIR)

#==============================================================================
clean: clean_objects clean_bin clean_lib clean_doc #clean_lattice

clean_objects:
	rm -rf $(OBJ_DIR)
	rm -f $(PRO_DIR)/*.o

clean_bin:
	rm -rf $(BIN_DIR)

clean_lib:
	rm -rf $(LIB_DIR)

clean_doc:
	rm -rf doc/html
	rm -rf doc/latex

clean_lattice:
	cd latticetester; make clean

mkdir:
	mkdir -p bin
	mkdir -p obj
	mkdir -p lib

#==============================================================================

.PHONY: clean latticetester

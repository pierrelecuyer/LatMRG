# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++14 -Wall -O2
DEBUG_FLAGS = -std=c++14 -g -Wall -O2

# The header files of the LatMRG library and the LatticeTester library
INCLUDES = -I./include -I./latticetester/include

#Definition to compile with yafu on if it is included
ifeq ($(wildcard data/yafu),)
    YAFU =
else
    YAFU = -DUSE_YAFU
endif

# Library path. This assumes NTL is in /usr/local/lib (its default path).
STAT_LIBS_PATH = -Wl,-Bstatic -L$(LIB_DIR)
STAT_LIBS = -llatmrg -llatticetester
DYN_LIBS_PATH = -Wl,-Bdynamic -L/usr/local/lib
DYN_LIBS = -lntl -lgmp -ltinyxml2 #-ltestu01 -lmylib

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
#MRG_HIGH  = $(SRC_DIR)/latmrg-high
#MRG_TYPES = $(SRC_DIR)/latmrg/mrgtypes

# The source files are in SRC_DIR. This grabs subdirectories
SRCS = $(wildcard $(MRG_LOW)/*.cc)
PROGS_CC = $(filter-out ./progs/SeekMain.cc, $(wildcard $(PRO_DIR)/*.cc))
PROGS_H = $(wildcard $(PRO_DIR)/*.h)
EX_CC = $(wildcard $(EX_DIR)/*.cc)

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
STA_OBJS = $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
DYN_OBJS = $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.so)
PROGS_O = $(PROGS_CC:$(PRO_DIR)/%.cc=$(PRO_DIR)/%.o)
EX_O = $(EX_CC:%.cc=%.o)

LATTEST_DEP = $(wildcard latticetester/progs/*.cc) \
	      $(wildcard latticetester/src/*.cc) \
	      $(wildcard latticetester/include/latticetester/*.h)


# A separator to segment the information printed on screen

SEP = @echo ================================================================================

default: $(PROGS_O)
	@echo
	@echo 'LatMRG programs compiled in ./bin'

all: clean_all default examples doc

#===============================================================================
# Building the library

$(LIB_DIR)/liblatmrg.a: $(STA_OBJS) | message_lib $(LIB_DIR)
	rm -f $(LIB_DIR)/liblatmrg.a
	ar rcs $(LIB_DIR)/liblatmrg.a $(STA_OBJS)
	$(CC) -o $(LIB_DIR)/liblatmrg.so -shared $(DYN_OBJS)
	@echo
	@echo 'LatMRG library archive created: ./lib/liblatmrg.a'
	@echo

message_lib:
	$(SEP)
	@echo 'Building library(ies)'
	@echo
	@make -q $(LIB_DIR)/liblatticetester.a && echo Everything up to date! || true


$(LIB_DIR):
	mkdir -p $(LIB_DIR)
	@echo

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/latmrg
	@echo

$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cc $(INC_DIR)/%.h | message_obj $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(YAFU) -c $< -o $@
	$(CC) -fPIC $(CFLAGS) $(INCLUDES) $(YAFU) -c $< -o $(@:%.o=%).so

message_obj:
	$(SEP)
	@echo 'Compiling the LatMRG library'
	@echo

#===============================================================================
# Building the programs of ./progs

$(PRO_DIR)/%.o: $(LIB_DIR)/liblatticetester.a $(LIB_DIR)/liblatmrg.a\
  $(PRO_DIR)/%.cc $(PROGS_H) | message_progs $(BIN_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(YAFU) -c $(PRO_DIR)/$(@:progs/%.o=%).cc -o $@
	$(CC) $@ $(STAT_LIBS_PATH) $(STAT_LIBS) $(DYN_LIBS_PATH) $(DYN_LIBS) \
	  -o $(BIN_DIR)/$(@:progs/%.o=%)

message_progs:
	$(SEP)
	@echo 'Building the LatMRG executables'

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
	@echo

#===============================================================================
# Building the documentation

doc:
	$(SEP)
	@echo 'Building the documentation'
	@echo
	doxygen docs/doc/doc-gen
	@echo 'Documentation now in ./docs'
	@echo

#===============================================================================
# Installation/removal of LatMRG

#===============================================================================
# Execution of the examples

PROGS_INPUTS = $(EX_DIR)/inputs

MK_XML = $(wildcard $(PROGS_INPUTS)/mk/mk*.xml)
PERIOD_XML = $(wildcard $(PROGS_INPUTS)/period/period*.xml)
LAT_XML = $(wildcard $(PROGS_INPUTS)/lat/lat*.xml)
SEEK_XML = $(wildcard $(PROGS_INPUTS)/seek/seek*.xml)

# Don't call this unless you are reckless
examples:default
	@echo
	$(SEP)
	@echo 'FindMK examples'
	./bin/MRGLattice $(MK_XML)
	@echo
	$(SEP)
	@echo 'Period examples'
	./bin/MRGLattice $(PERIOD_XML)
	@echo
	$(SEP)
	@echo 'Lattest examples'
	./bin/MRGLattice $(LAT_XML)
	@echo
	$(SEP)
	@echo 'Seek examples'
	./bin/MRGLattice $(SEEK_XML)

clean_ex:
	rm -f $(PROGS_INPUTS)/*/*.res

#===============================================================================
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

$(LIB_DIR)/liblatticetester.a:$(LATTEST_DEP) | message_lib $(LIB_DIR)
	$(SEP)
	@echo 'Building/Updating LatticeTester with waf'
	@echo
	cd latticetester; ./waf build
	@echo 'LatticeTester build/update finished'
	@echo
	cp latticetester/build/src/liblatticetester.a $(LIB_DIR)
	@echo

clean_lattice:
	$(SEP)
	@echo 'Cleaning LatticeTester'
	@echo
	cd latticetester; ./waf clean
	@echo

#===============================================================================
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
	rm -rf docs/*.html
	rm -rf docs/*.js
	rm -rf docs/*.png
	rm -rf docs/*.css
	rm -rf docs/search

#===============================================================================
# Graphical stuff

separator:
	$(SEP)
	@echo

#===============================================================================
# PHONY targets

.PHONY: doc clean clean_all examples $(EX_BUILD)/ $(EX_CC) separator\
  config_latticetester $(MK_DAT) mk_ex_head seek_ex_head all_ex

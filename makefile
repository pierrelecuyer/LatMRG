#
# To use this makefile, you have to:
#       - [line 19, 24] Modifie include and library paths
#       - [line 32} Modifie the path for main.cc if it is not in src folder, and update .cc list
# 'make depend' uses makedepend to automatically generate dependencies
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++11 -g -Wall

# define any directories containing header files other than /usr/include
##INCLUDES = -I../include -I/u/simul/opt/boost-1.60.0/include -I/u/jemelaym/Code/ntl/include
INCLUDES = \
-I/Users/Erwan1/projects/github/LatMRG/include \
-I/Users/Erwan1/projects/github/Lattice\ Tester/include \
-I/Users/Erwan1/projects/github/Lattice\ Tester \
-I/usr/local/include \
-I/usr/local/Cellar/NTL \
-I/usr/local/Cellar/boost/1.60.0_2/include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = \
-L/usr/local/lib/ \
-L/usr/local/Cellar/boost/1.60.0_2/lib/

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lntl -lgmp -lm -ltestu01 -lmylib

# define the C source files
SRCS = LatMain.cc \
\
../Lattice\ Tester/SimpleMRG.cc \
../Lattice\ Tester/src/IntLatticeBasis.cc \
../Lattice\ Tester/src/Const.cc \
../Lattice\ Tester/src/CoordinateSets.cc \
../Lattice\ Tester/src/Coordinates.cc \
../Lattice\ Tester/src/IntFactor.cc \
../Lattice\ Tester/src/Lacunary.cc \
../Lattice\ Tester/src/NormaBestLat.cc \
../Lattice\ Tester/src/NormaLaminated.cc \
../Lattice\ Tester/src/NormaMinkL1.cc \
../Lattice\ Tester/src/NormaMinkowski.cc \
../Lattice\ Tester/src/NormaPalpha.cc \
../Lattice\ Tester/src/NormaRogers.cc \
../Lattice\ Tester/src/Normalizer.cc \
../Lattice\ Tester/src/Num.cc \
../Lattice\ Tester/src/OrderDependentWeights.cc \
../Lattice\ Tester/src/PODWeights.cc \
../Lattice\ Tester/src/ProductWeights.cc \
../Lattice\ Tester/src/ProjectionDependentWeights.cc \
../Lattice\ Tester/src/Random.cc \
../Lattice\ Tester/src/Reducer.cc \
../Lattice\ Tester/src/Util.cc \
\
./src/AverageCaseMerit.cc \
./src/HyperplaneDistanceMerit.cc \
./src/LatticeTest.cc \
./src/PalphaLCG.cc \
./src/ProjIteratorSuccCoords.cc \
./src/SpectralMeritClassic.cc \
./src/BoundJSInter.cc \
./src/IntFactorization.cc \
./src/MRGComponent.cc \
./src/PalphaOrderDependent.cc \
./src/Rank1Lattice.cc \
./src/Table.cc \
./src/BoundJSstar.cc \
./src/IntLattice.cc \
./src/MRGComponentFactory.cc \
./src/PalphaProduct.cc \
./src/ReportFooterLat.cc \
./src/TestProjections.cc \
./src/Chrono.cc \
./src/IntPrimitivity.cc \
./src/MRGLattice.cc \
./src/ParamReader.cc \
./src/ReportHeaderLat.cc \
./src/WorstCaseMerit.cc \
./src/CoefApproxFact.cc \
./src/KorobovLattice.cc \
./src/MRGLatticeFactory.cc \
./src/ParamReaderLat.cc \
./src/ReportLat.cc \
./src/Writer.cc \
./src/CoefEqual.cc \
./src/LatConfig.cc \
./src/MRGLatticeLac.cc \
./src/ParamReaderSeek.cc \
./src/Searcher.cc \
./src/WriterRes.cc \
./src/CoefZero.cc \
./src/LatTestAll.cc \
./src/Merit.cc \
./src/PolyPE.cc \
./src/SearcherCBC.cc \
./src/Zone.cc \
./src/Discrepancy.cc \
./src/LatTestBeyer.cc \
./src/Modulus.cc \
./src/Primes.cc \
./src/SearcherKorobov.cc \
./src/DoubleFormatter.cc \
./src/LatTestPalpha.cc \
./src/OrdepBound.cc \
./src/ProjIteratorDefault.cc \
./src/SeekConfig.cc \
./src/ExactDiscStar.cc \
./src/LatTestSpectral.cc \
./src/Palpha.cc \
./src/ProjIteratorNonSuccCoords.cc \
./src/ShortestDualVectorMerit.cc

#../latticetester/src/Basis.cc \
../latticetester/src/IntLattice.cc \
../latticetester/src/Rank1Lattice.cc \

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = main

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY:	depend clean

all:	$(MAIN)
	@echo  Simple compiler named mycc has been compiled

$(MAIN): $(OBJS)
	$(CXX) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
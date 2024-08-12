#ifndef LATMRG_EXECCOMMON_H
#define LATMRG_EXECCOMMON_H

/**
 * This file contains the header inclusions needed by the executables and defines
 * global variables with the default values for the program.
 *
 * */
#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Reducer.h"
#include "latticetester/Util.h"

#include "../include/latmrg/MWCComponent.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/ComboLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/Projections.h"
#include "latmrg/FigureOfMerit.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/Projections.h"
#include "latmrg/Test.h"

#ifndef MRGLATTICE_MAIN_EXEC
typedef NTL::ZZ Int;
typedef double Real;
typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;
typedef NTL::vector<Real> RealVec;
#endif

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace LatMRG;
using namespace LatticeTester::Random;

// Variable definitions for executables only
const std::string delim = "\n========================================"
"========================================\n\n";

#ifdef MRGLATTICE_MAIN_EXEC
// Types to use for computations
std::string types = "ZD";
#endif


#define LATMRG_SEEK
/**
 * This stores the configuration of a problem. This contains many parameters
 * used in the executables and when using the `test()` function.
 *
 * This structure is intended to be included in executables only and needs to
 * be compiled on demand. To do so, simply `#define LATMRG_TEST` before
 * including this header.
 * */
template<typename Integ, typename Real> struct ConfigSeek {
#include "Config.h"
};
/**
 * This performs a test on a lattice. That is, given a lattice and a `Config`
 * it populates a `FigureOfMerit` objects and returns it.
 * */
template<typename Lat>
FigureOfMerit<Lat> test_seek(Lat & lattice, ConfigSeek<typename Lat::Int, typename Lat::Real>& conf) {
#include "Test.h"
}
#undef LATMRG_SEEK
#define LATMRG_LAT
/**
 * This stores the configuration of a problem. This contains many parameters
 * used in the executables and when using the `test()` function.
 *
 * This structure is intended to be included in executables only and needs to
 * be compiled on demand. To do so, simply `#define LATMRG_TEST` before
 * including this header.
 * */
template<typename Integ, typename Real> struct ConfigLat {
#include "Config.h"
};
/**
 * This performs a test on a lattice. That is, given a lattice and a `Config`
 * it populates a `FigureOfMerit` objects and returns it.
 * */
template<typename Lat>
FigureOfMerit<Lat> test_lat(Lat & lattice, ConfigLat<typename Lat::Int, typename Lat::Real>& conf) {
#include "Test.h"
}
#undef LATMRG_LAT

#endif


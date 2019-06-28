#ifndef LATMRG_EXEC_H
#define LATMRG_EXEC_H

/**
 * This file contains functions, classes and code used in SeekRe. This file is
 * used to remove some of the clutter in SeekRe.cc. This is because some of the
 * functions needed already exist but need to be rewritten.
 * */
#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/Const.h"
#include "latticetester/Reducer.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/ComboLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/Projections.h"
#include "latmrg/FigureOfMerit.h"
#include "latmrg/Const.h"
#include "latmrg/Projections.h"

typedef NTL::ZZ Int;
typedef double Dbl;
typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;
typedef NTL::vector<Dbl> DblVec;

#define LATMRG_USE_CONFIG
#define LATMRG_USE_TEST
#include "latmrg/Test.h"

using LatticeTester::IntLattice;
using LatticeTester::Normalizer;

// Variable definitions for executables only
const std::string delim = "\n========================================"
"========================================\n\n";
#ifdef MRGLATTICE_MAIN_EXEC
LatMRG::GenType gen_type = LatMRG::MRG;
LatticeTester::CriterionType crit_type = LatticeTester::SPECTRAL;
bool dual = true;
LatticeTester::PreReductionType red_type = LatticeTester::FULL;
LatticeTester::NormaType norma_type = LatticeTester::NONE;
double time_limit = 3600.0;
int max_gen = 50;
int detail = 1;
bool period = false;

// Values without default
LatMRG::Projections* proj;
std::string gen_string;

#endif


#endif

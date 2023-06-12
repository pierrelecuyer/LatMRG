#ifndef LATMRG_CONFIGSEEK_H
#define LATMRG_CONFIGSEEK_H

#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/Reducer.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/ComboLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/Projections.h"
#include "latmrg/FigureOfMerit.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/Projections.h"
#include "latmrg/Test.h"

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace LatMRG;
using namespace LatticeTester::Random;

typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

template<typename Int, typename Real> struct ConfigSeek {


	#include "Config.h"
};







// Type of figure of merit
LatticeTester::NormaType normaType = LatticeTester::NONE;
LatticeTester::CriterionType criterion = LatticeTester::SPECTRAL;
LatticeTester::PreReductionType reduction = LatticeTester::FULL;
LatticeTester::NormType norm = LatticeTester::L2NORM;



bool use_dual = true;    // Must be clearly defined!
#ifdef LATMRG_SEEK
bool best = true;
int num_gen = 0;
Real currentMerit = Real(0);
std::string construction = "RANDOM";
#endif
// Projection
Projections* proj;



// MRGPeriod stores the stuff we might want to know, such as the modulo
// the order and even the coefficients
int num_comp = 1;
std::vector<std::string> search_mode;
std::vector<bool> period;
std::vector<MRGPeriod<Int>*> fact;
// This is used to store the coefficients of a MRG and the information on how to
// search for coefficients of a generator. This can have multiple different formats.
std::vector<IntVec> coeff;

bool gen_set = false;
bool test_set = false;
bool proj_set = false;

// Program variables
double timeLimit = 600;


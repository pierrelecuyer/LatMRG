/**
 * In this example, we show that the results of the new LatMRG version
 * are consistent with the examples reported in the the article [rLEC97c].
 * It is important that there are either no difference or, if they appear,
 * it needs to be shown that the calculations in the current version
 * are correct.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <algorithm>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/ReducerStatic.h"

#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGLatticeLac.h"
#include "latmrg/MLCGLattice.h"
#include "latmrg/MLCGLatticeLac.h"

using namespace LatticeTester;

/**
 * This function calculates the FOM of the dual lattice of an MRG with modulus 'm'
 * and vector of multipliers 'aa', for the successive coordinates in `lowDim` to `highDim` dimensions.
 */
template<typename Int, typename Real>
static void FOMSuccLattice(Int &m, IntVec &aa, int64_t lowDim, int64_t highDim) {
   int64_t k = aa.length() - 1;  // order
   NTL::Vec<int64_t> t;
   LatMRG::MRGLattice<Int, Real> mrg(m, aa, highDim);
   ReducerBB<Int, Real> red(mrg);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLattice: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);
   if (k == 1) {
      LatMRG::LCGLattice<Int, Real> lcg(m, aa[1], highDim);
      red.setIntLattice(lcg);
      std::cout << "\nUsing LCGLattice: \n";
      std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
      fomdual.computeMeritSucc(lcg, lowDim, highDim);
   }
}

/**
 * Same as the previous function, but with lacunary indices defined by 'lac'.
 */
template<typename Int, typename Real>
static void FOMSuccLatticeLac(Int &m, IntVec &aa, IntVec lac, int64_t lowDim, int64_t highDim) {
   int64_t k = aa.length() - 1;  // order
   NTL::Vec<int64_t> t;
   LatMRG::MRGLatticeLac<Int, Real> mrg(m, aa, highDim);
   mrg.setLac(lac);
   ReducerBB<Int, Real> red(mrg);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLatticeLac: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);
   if (k == 1) {
    LatMRG::LCGLatticeLac<Int, Real> lcg(m, aa[1], highDim);
    lcg.setLac(lac);
    red.setIntLattice(lcg);
    std::cout << "\nUsing LCGLatticeLac: \n";
    std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
    fomdual.computeMeritSucc(lcg, lowDim, highDim);
    }
}

int main() {

   int64_t maxdim = 12;
   std::cout << "Types: NTL::ZZ, double \n";
   // NTL::Vec<int64_t> t; // The t-vector for the FoM.

   Int m;
   IntVec aa; // Vector a has size k+1, if it contains k elements.

   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 2147483647 and a = 45991\n";
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 45991;
   FOMSuccLattice<Int, Real>(m, aa, 2, maxdim);
   std::cout << "\nThe results should be compared with Table 6 of [rLEC97c].\n";

   std::cout << "\n=============================================================\n";
   std::cout << "An MRG with m = 9223372036854773561 and"
         << " aa = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   aa.SetLength(4);
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   FOMSuccLattice<Int, Real>(m, aa, 4, maxdim);
   std::cout << "\nThe results should be compared with Table 7 of [rLEC97c].\n";

   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 4294967296, a = 1099087573. \n";
   maxdim = 30;
   m = 4294967296;
   m = m / 4;
   aa.SetLength(2);
   aa[1] = 1099087573;
   FOMSuccLattice<Int, Real>(m, aa, 2, maxdim);
   std::cout << "\nThe results and timings should be compared with Table 1 of [rLEC97c].\n";

   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 2147483647, a = 16807, lacunary with leap d = 131072, \n";
   std::cout << "Lacunary indices 1, 2, 3, d+1, d+2, d+3, 2*d+1, 2*d+2, ...\n";
   maxdim = 42;
   Int leap(131072);  // Value of d
   IntVec lac;
   lac.SetLength(maxdim);
   for (int64_t i = 0; i < maxdim; i++)
      lac[i] = 1 + (i % 3) + leap * (i / 3);
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 16807;
   FOMSuccLatticeLac<Int, Real>(m, aa, lac, 2, maxdim);
   std::cout << "\nThe results should be compared with Table 2 of [rLEC97c].\n";
}

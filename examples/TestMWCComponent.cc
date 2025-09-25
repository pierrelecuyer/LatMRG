/**
 * TestMWCComponent: Testing MWC generators.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/Random.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"

using namespace LatticeTester;
using namespace LatMRG;

template<typename Int>
static void tryMultipliers(int64_t k, int64_t e, int64_t alpha, int64_t numMult, int64_t maxdim) {
   Int b = getPow2<Int>(e);
   std::cout << "Order k = " << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   std::cout << "We try " << numMult << " random vectors of coefficient with " << alpha << " bits each.\n";
   std::cout << "We found the following candidates: \n\n";
   Int m;
   IntVec aa;
   aa.SetLength(k + 1);
   aa[0] = Int(1);
   LCGLattice<Int, Real> lcg(maxdim);
   ReducerBB<Int, Real> red(lcg);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(maxdim);
   NTL::Vec<int64_t> t;
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   for (int64_t i = 0; i < numMult; i++) {
      for (int64_t j = 1; j <= k; j++)
         aa[j] = conv<Int>(LatticeTester::RandBits(alpha));
      m = computeLCGModulusMWC(b, aa);
      if (mIsPrime(m) <= 1) {  // Prime or probably prime.
         // std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
         // std::cout << "  vector aa = " << aa << "\n";
         if (mIsSafePrime(m, 50) <= 1) {
            std::cout << "i = " << i << ", Modulo m = " << m << "  is probably a safe prime.\n";
            IntFactorization<Int> fact(m - Int(1));
            LatMRG::PrimeType ptype = fact.decompToFactorsInv(DECOMP, NULL);
            // std::cout << "  decomp done, num factors: " <<  (m_fact.getFactorList()).size() << "\n";
            if (ptype > 1) {
               std::cout << "  Factorization failed. " << "\n";
            }
            if (isPrimitiveElement(Int(2), fact, m - Int(1)))
               std::cout << "  We also find that 2 is a primitive element mod m.  Spectral test:\n";
            // Since `m` has changed, we must recompute the normalization bounds.
            lcg.setModulus(m);
            lcg.seta(b, maxdim);
            normaDual.computeBounds(-log(m), 1);
            fomdual.computeMeritSucc(lcg, k+1, maxdim);
            std::cout << "\n";
         }
      }
   }
   tmp = clock() - tmp;
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: ";
   std::cout << "We make a search for MWC generators of order k, modulo b, a_0 = 1,\n";
   std::cout << "period length (m-1)/2, and we apply the spectral test (in dual) to the retained ones.\n\n";

   int64_t k = 3;
   int64_t e = 64;  // b = 2^{64}
   int64_t alpha = 62;   // Number of bits for the random coefficients.
   int64_t maxdim = 12;  // Max dimension for lattice.
   tryMultipliers<Int>(k, e, alpha, 100000, maxdim);
   return 0;
}

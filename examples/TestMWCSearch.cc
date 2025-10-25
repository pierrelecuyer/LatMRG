/**
 * TestMWCSearch: Search and compare MWC generators with different constraints.
 * This is an experimental program, incomplete ....
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
#include "latmrg/MWCLattice.h"

using namespace LatticeTester;
using namespace LatMRG;

int64_t maxdim = 8;
LCGLattice<Int, Real> lcg(maxdim);
MWCLattice<Int, Real> mwc(maxdim);
ReducerBB<Int, Real> red(mwc);
WeightsUniform weights(1.0);
NormaBestLat normaDual(maxdim);
FigureOfMeritDualM<Int, Real> fomdual(weights, normaDual, &red);
double bestmerit = 0.0;   // Best merit value so far.


template<typename Int, typename Real>
static double testSafePrime(Int m, Int b, IntVec aa) {
   int64_t k = aa.length()-1;
   std::cout << "Coefficients aa = " << aa << "\n";
   std::cout << " in hexadecimal = ";
   for (int64_t j = 0; j <= k; j++)
      std::cout << " " << "0x" << std::hex << conv<uint64_t>(aa[j]);
   std::cout << std::dec << "\n";
   std::cout << "b - a_k = " << b - aa[k] << "\n";
   Int sum = b;
   for (int64_t j = 1; j <= k; j++)
       sum -= aa[j];
   std::cout << "b - sum a_j = " << b - sum << "\n";
   std::cout << "Modulo m = " << m << "  is a safe prime.\n";
   std::cout << "Log_2(m) = " << Lg(m) << "\n";

   IntFactorization<Int> fact(m - Int(1));
   LatMRG::PrimeType ptype = fact.decompToFactorsInv(DECOMP, NULL);
   if (ptype > 1) std::cout << "  Factorization failed. " << "\n";
   if (isPrimitiveElement(Int(2), fact, m))
      std::cout << "We also find that 2 is a primitive element mod m. \nSpectral test:\n";
   Int suma2 = aa[0] * aa[0];
   for (int64_t j=1; j <= k; j++)  suma2 += aa[j] * aa[j];
   std::cout << "a_0^2 + ... + a_k^2 = " << conv<double>(suma2) << "\n";
   std::cout << "bound 2 dim: (b^2 + 1)^{-1/2} = " << 1.0 / sqrt(conv<double>(b*b + 1)) << "\n";
   std::cout << "bound k+1 dim: (a_0^2 + ... + a_k^2)^{-1/2} = " << 1.0 / sqrt(conv<double>(suma2)) << "\n";
   std::cout << "spectral test in dual lattice:\n";
   // lcg.setModulus(m);
   // lcg.seta(b, maxdim);
   normaDual.computeBounds(-log(m), 1);
   // fomdual.computeMeritSucc(lcg, k+1, maxdim);  // Start in k+1 dim.

   mwc.setb(b);
   mwc.setaa(aa, maxdim);
   double merit = fomdual.computeMeritSucc(mwc, k+1, maxdim);  // Start in k+1 dim.
   //mwc.buildDualBasis(k+2);
   //std::cout << " Initial dual basis:\n " << mwc.getDualBasis() << "\n";
   std::cout << "\n";
   return merit;
}

// Order k, b=2^e, er random bits for the a_j.
template<typename Int, typename Real>
static void tryMultipliers(int64_t k, int64_t e, int64_t er, int64_t numMult) {
   double merit;
   double bestmerit = 0.0;
   IntVec bestaa;
   bestaa.SetLength(k + 1);

   Int b = NTL::power(Int(2), e);
   // Int sqrtb = NTL::power(Int(2), e/2);
   std::cout << "Order k = " << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   std::cout << "We want b - a_k < 2^{er}, with er = " << er << "\n";
   std::cout << "We try " << numMult << " random vectors of coefficient.\n";
   std::cout << "We found the following candidates: \n\n";
   Int m;
   IntVec aa;
   aa.SetLength(k + 1);
   // aa[0] = -1;
   aa[0] = -Int(1) - NTL::power(Int(2), 31);
   // aa[0] = conv<Int>(LatticeTester::RandBits(60));
   for (int64_t j = 1; j < k-1; j++) aa[j] = 0;
   // fomdual.setVerbosity(4);
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   for (int64_t i = 1; i <= numMult; i++) {
      // aa[k] = b - conv<Int>(i);
      aa[k] = conv<Int>(LatticeTester::RandBits(62));
      //aa[k] = b - NTL::power(Int(2), 60) + conv<Int>(LatticeTester::RandBits(58));
      // aa[k-1] = conv<Int>(LatticeTester::RandBits(58));
      // aa[k-2] = conv<Int>(LatticeTester::RandBits(58));
      // aa[k-2] = b - NTL::power(Int(2), er) + conv<Int>(LatticeTester::RandBits(er));
      // aa[k] = b/3 + conv<Int>(LatticeTester::RandBits(e-1));
         // aa[j] = b - sqrtb + conv<Int>(LatticeTester::RandBits(e/2-1));
         //aa[j] = b/2 + conv<Int>(LatticeTester::RandBits(e-1));
         // aa[j] = (1 << (alpha-1)) + conv<Int>(LatticeTester::RandBits(alpha-1));
      // aa[k-1] = aa[k];
      // aa[k-2] = aa[k-1];
      m = computeLCGModulusMWC(b, aa);
      if (mIsPrime(m) <= 1) {  // Prime or probably prime.
         // std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
         // std::cout << "  vector aa = " << aa << "\n";
         if (mIsSafePrime(m, 50) <= 1) {
            std::cout << "i = " << i << "\n";
            merit = testSafePrime<Int, Real>(m, b, aa);
            if (merit > bestmerit) {
               bestmerit = merit;
               bestaa = aa;
            }
         }
      }
   }
   std::cout << "\n\n";
   std::cout << "----------------------------------------------------------------\n";
   std::cout << "Best FOM: " << bestmerit << "\n";
   std::cout << "Best aa: " << bestaa << "\n\n";
   m = computeLCGModulusMWC(b, bestaa);
   merit = testSafePrime<Int, Real>(m, b, bestaa);
   std::cout << "Merit = " << merit << "\n";
   std::cout << "----------------------------------------------------------------\n\n";

   tmp = clock() - tmp;
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
}


int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: ";
   std::cout << "We make a search for MWC generators of order k, modulo b, \n";
   std::cout
         << "period length (m-1)/2, and we apply the spectral test (in dual) to the retained ones.\n\n";
   fomdual.setVerbosity(4, 0.01);

   int64_t k = 3;
   int64_t e = 64;  // b = 2^{64}
   int64_t er = 64;   // Number of random bits for the random coefficients.
   // int64_t maxdim = 12;  // Max dimension for lattice.
   tryMultipliers<Int, Real>(k, e, er, 1000000);

   return 0;

}

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

int64_t maxdim = 12;
int64_t maxdimproj = 5;
LCGLattice<Int, Real> lcg(maxdim);
MWCLattice<Int, Real> mwc(maxdim);
IntLattice<Int, Real> proj(maxdimproj);  // For the projections.
ReducerBB<Int, Real> red(mwc);
WeightsUniform weights(1.0);
NormaBestLat normaDual(maxdim);
FigureOfMeritDualM<Int, Real> fomdual(weights, normaDual, &red);

double bestmerit = 0.0;   // Best merit value so far.
double merit;
int64_t numSafe = 0;
// int64_t numOrderp = 0;
int64_t numPrimitive = 0;
int64_t verbose = 0;

// Tests a given set of parameters for which `m` is a safe prime,
// for successive values in up to `maxdim` dimensions, and prints info on it.
// If `nonsucc`, also tests for non-successive projections.
template<typename Int, typename Real>
static double testSafePrime(Int m, Int b, IntVec aa, bool nonsucc) {
   int64_t k = aa.length() - 1;
   Int a0inv = InvMod(-aa[0], b);
   if (((-aa[0] * a0inv) % b) != 1) std::cout << "Error: This inverse is wrong! \n";
   // std::cout << std::dec << "\n";
   Int sum = b;
   for (int64_t j = 1; j <= k; j++)
      sum -= aa[j];
   if (verbose > 0) {
      std::cout << "Coefficients aa = " << aa << "\n";
      std::cout << " in hexadecimal = ";
      std::cout << " " << "- 0x" << std::hex << conv<uint64_t>(-aa[0]);
      for (int64_t j = 1; j <= k; j++)
         std::cout << " " << "0x" << std::hex << conv<uint64_t>(aa[j]);
      cout << "\nInverse of " << -aa[0] << " mod b is " << a0inv << " = 0x" << std::hex
            << conv<uint64_t>(a0inv) << endl;
      std::cout << "b - a_k = " << b - aa[k] << "\n";
      std::cout << "b - sum a_j = " << b - sum << "\n";
      std::cout << "Modulo m = " << m << "  is a safe prime.\n";
      std::cout << "Log_2(m) = " << Lg(m) << "\n";
      Int suma2 = aa[0] * aa[0];
      for (int64_t j = 1; j <= k; j++)
         suma2 += aa[j] * aa[j];
      std::cout << "Spectral test:\n";
      std::cout << "a_0^2 + ... + a_k^2 = " << conv<double>(suma2) << "\n";
      std::cout << "bound 2 dim: (b^2 + 1)^{-1/2} = " << 1.0 / sqrt(conv<double>(b * b + 1))
            << "\n";
      std::cout << "bound k+1 dim: (a_0^2 + ... + a_k^2)^{-1/2} = "
            << 1.0 / sqrt(conv<double>(suma2)) << "\n";
      std::cout << "spectral test in dual lattice:\n\n" << std::dec;
   }
   lcg.setModulus(m);
   lcg.seta(b, maxdim);
   normaDual.computeBounds(-log(m), 1);
   // fomdual.computeMeritSucc(lcg, k+1, maxdim);  // Start in k+1 dim.
   //mwc.setb(b);
   //mwc.setaa(aa, maxdim);   // must implement buildProjection  ****************
   double merit = fomdual.computeMeritSucc(lcg, k + 1, maxdim);  // Start in k+1 dim.
   // std::cout << " min merit succ = \n " << merit << "\n";
   //mwc.buildDualBasis(k+2);
   // std::cout << " Initial dual basis:\n " << mwc.getDualBasis() << "\n";
   if (nonsucc) {
      if (verbose > 0) std::cout << "\n";
      proj.setModulus(m);
      merit = min(merit, fomdual.computeMeritNonSucc(lcg, proj));
   }
   return merit;
}

// Tests if m is a safe prime and has max order if a_0 < -1.
template<typename Int>
static bool testOrder (Int &m, Int &b, IntVec &aa) {
   if (mIsPrime(m) > 1) return false;
   // Now we know m is prime or probably prime.
   if (IntFactor<Int>::isPrime((m - Int(1)) / Int(2), 50)) {
      // We have a safe prime.
      numSafe++;
      if (aa[0] == -1) return true;
      if (((m % Int(8)) == 3) | ((m % Int(8)) == 5)) {
         numPrimitive++;
         return true;
      }
   }
   return false;
}

// Tests if 2 and b are primitive elements mod m.
template<typename Int, typename Real>
static bool testPrimitive (Int &m, Int &b, IntVec &aa) {
   IntFactorization<Int> fact(m - Int(1));
   LatMRG::PrimeType ptype = fact.decompToFactorsInv(DECOMP, NULL);
   if (ptype > 1) std::cout << "  Factorization failed. " << "\n";
   bool primb = isPrimitiveElement(b, fact, m);
   bool prim2 = isPrimitiveElement(Int(2), fact, m);
   //if (verbose > 0)
   std::cout << "Is 2 primitive mod m? " << prim2 << "\n";
   //if (verbose > 0)
   std::cout << "Is b primitive mod m? " << primb << "\n";
   return primb;
}


/**
 * Order k, b=2^e, `numaj` nonzero coefficients `a_j` besides a_0,
 * e0 random bits for a_0, ek random bits for a_k, ej random bits for other coefficients,
 * if `nonsucc` we look at non-succ projections in the FOM.
 * We try `numMult` random vectors of parameters for a MWC with given constraints.
 * We test all those for which `m` is a safe prime, and finally we print the best FOM.
 */
template<typename Int, typename Real>
static void tryMultipliers(int64_t k, int64_t numaj, int64_t e, int64_t e0, int64_t ek, int64_t ej,
      int64_t numMult, bool nonsucc) {
   assert(numaj <= k);
   numSafe = 0;
   // numOrderp = 0;
   numPrimitive = 0;
   bestmerit = 0.0;
   IntVec bestaa;
   bestaa.SetLength(k + 1);
   Int b = NTL::power(Int(2), e);
   Int m;
   IntVec aa;
   aa.SetLength(k + 1);
   aa[0] = -1;
   for (int64_t j = 1; j < k - 1; j++)
      aa[j] = 0;

   // std::cout << "e = " << e << "\n";
   // Int sqrtb = NTL::power(Int(2), e/2);
   std::cout << "Order k = " << std::dec << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   if (e0 == 0) std::cout << "We impose a_0 = -1." << "\n";
   std::cout << "Number of nonzero coefficients a_j besides a_0: " << numaj << "\n";
   std::cout << "Random bits for coefficients, e_0 = " << std::dec << e0 << ",  "
         "e_k = " << ek << ", e_j = " << ej << "\n";
   std::cout << "We try " << numMult << " random vectors of coefficient.\n";
   if (verbose > 0) std::cout << "We found the following candidates: \n\n";
   // fomdual.setVerbosity(4);
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   for (int64_t i = 1; i <= numMult; i++) {
      if (e0 > 0) {
         aa[0] = -conv<Int>(LatticeTester::RandBits(e0));
         if ((aa[0] % 2) == 0) aa[0] += 1;
      }
      // aa[k] = b - conv<Int>(i);   // If we want a_k very close to b, for testing.
      // aa[k] = b - NTL::power(Int(2), 60) + conv<Int>(LatticeTester::RandBits(58));
      aa[k] = conv<Int>(LatticeTester::RandBits(ek));
      if (numaj > 1) aa[k - 1] = conv<Int>(LatticeTester::RandBits(ej));
      if (numaj > 2) aa[k - 2] = conv<Int>(LatticeTester::RandBits(ej));
      // aa[k-1] = aa[k];  // If we want equal coefficients, for testing.
      // aa[k-2] = aa[k-1];
      m = computeLCGModulusMWC(b, aa);
      if (m < Int(2)) std::cout << "tryMultipliers: Modulo m too small, m = " << m << "\n";
      if (testOrder<Int>(m, b, aa)) {
         if (verbose > 0) std::cout << "i = " << std::dec << i << ", testing a safe prime \n";
         merit = testSafePrime<Int, Real>(m, b, aa, false);
         if (merit > bestmerit) {
             bestmerit = merit;
             bestaa = aa;
         }
      }
   }
   if (bestmerit == 0.0) {
      std::cout << "No candidate with `m` safe prime was found. \n\n\n";
      return;
   }
   // Prints info on the best candidate.
   if (verbose > 0) std::cout << "\n\n";
   std::cout << "----------------------------------------------------------------\n";
   std::cout << "Number safe primes: " << numSafe << "\n";
   if (aa[0] < -1) std::cout << "Number primitive: " << numPrimitive << "\n";
   std::cout << "Best FOM: " << bestmerit << "\n";
   std::cout << "Best aa: " << bestaa << "\n\n";
   m = computeLCGModulusMWC(b, bestaa);
   verbose = 2;
   fomdual.setVerbosity(4, 0.01);
   merit = testSafePrime<Int, Real>(m, b, bestaa, nonsucc);
   verbose = 0;
   fomdual.setVerbosity(2, 0.01);
   std::cout << "Merit = " << merit << "\n";
   tmp = clock() - tmp;
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
   std::cout << "----------------------------------------------------------------\n\n";
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: ";
   std::cout << "We make a search for MWC generators of order k, modulo b, \n";
   std::cout << "period (m-1)/2, and apply spectral test (in dual) to retained ones.\n\n";
   fomdual.setVerbosity(2, 0.01);   // Do not print the details for all safe primes.
   verbose = 0;

   NTL::Vec<int64_t> t; // The t-vector for the FOM.
   t.SetLength(5);
   t[0] = 12;   // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 0;    // pairs
   t[2] = 8;    // triples
   t[3] = 6;    // etc.
   t[4] = 5;
   fomdual.setTVector(t, true);

   int64_t k = 1;
   int64_t numaj = 1;
   int64_t e = 64;    // b = 2^{64}
   int64_t e0 = 60;   // Number of random bits for a_0. If 0, then a_0 = -1.
   int64_t ek = 62;   // Number of random bits for a_k.
   int64_t ej = 58;   // Number of random bits for other a_j's.
   // int64_t maxdim = 12;  // Max dimension for lattice.
   int64_t numMult = 1000 * 10000;  // 10 millions

   k = 1;
   // tryMultipliers<Int, Real>(1, numaj, 6, 5, 5, 4, 25, true);
   tryMultipliers<Int, Real>(k, numaj, e, e0, ek, ej, numMult, true);

   k = 2;
   t[1] = 0;  // No pairs
   fomdual.setTVector(t, true);
   tryMultipliers<Int, Real>(2, 1, 64, 0, 62, 58, numMult, true);   // MWC1k2-0-62.res
   tryMultipliers<Int, Real>(2, 2, 64, 0, 62, 58, numMult, true);   // MWC1k2-62-60.res

   tryMultipliers<Int, Real>(2, 1, 64, 60, 62, 58, numMult, true);  // MWCa0k2-0-62.res
   tryMultipliers<Int, Real>(2, 2, 64, 60, 62, 58, numMult, true);  // MWCa0k2-62-60.res

   k = 3;
   t[2] = 0;  // No triples
   fomdual.setTVector(t, true);
   tryMultipliers<Int, Real>(3, 1, 64, 0, 62, 58, numMult, true);   // MWC1k3-00-62.res
   tryMultipliers<Int, Real>(3, 2, 64, 0, 60, 58, numMult, true);   // MWC1k3-0-60-58.res
   tryMultipliers<Int, Real>(3, 3, 64, 0, 60, 58, numMult, true);   // MWC1k3-60-58.res

   tryMultipliers<Int, Real>(3, 1, 64, 60, 62, 58, numMult, true);  // MWCa0k3-00-62.res
   tryMultipliers<Int, Real>(3, 2, 64, 60, 60, 58, numMult, true);  // MWCa0k3-0-60-58.res
   tryMultipliers<Int, Real>(3, 3, 64, 60, 60, 58, numMult, true);  // MWCa0k3-60-58.res

   return 0;

}

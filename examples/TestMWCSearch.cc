/**
 * TestMWCSearch: For various choices of `k` and `b`, and certain constraints on
 * the parameters, this program samples vectors of coefficients `aa` at random
 * and retains the MWC with best FOM among those having the maximal period `(m-1)/2`.
 * The FOM is based on the spectral test in the m-dual lattice.  The shortest vectors
 * are computed exactly using BB after a pre-reduction with BKZ with default parameters.
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

int64_t maxdim = 12;     // The maximal dimension.
int64_t maxdimproj = 5;
LCGLattice<Int, Real> lcg(maxdim);
MWCLattice<Int, Real> mwc(maxdim);
IntLattice<Int, Real> proj(maxdimproj);  // For the projections.
ReducerBB<Int, Real> red(mwc);
WeightsUniform weights(1.0);
NormaBestLat normaDual(maxdim);
FigureOfMeritDualM<Int, Real> fomdual(weights, normaDual, &red);
NTL::Vec<int64_t> t; // The t-vector for the FOM.

double bestmerit = 0.0;   // Best merit value so far.
double merit;             // For current candidate.
bool nonsucc = false;     // For projections onto non-successive coordinates.
bool onlySafe = true;     // If `true`, we retain only the safe primes,
int64_t numSafe = 0;      // Number of safe primes found.
int64_t numMaxPer = 0;    // Number of cases with period p = (m-1)/2.
int64_t verbose = 0;      // Selects the level of verbosity.

/**
 * Prints information and applies the spectral test for a given MWC instance
 * with given MWC modulus `b`, vectors of coefficients is `aa`, and LCG modulus `m`.
 * The spectral test is applied for successive coordinates in dimensions
 * from `mindim` up to `maxdim`. If `nonsucc` is true, it is also applied for
 * non-successive projections according to the vector `t`.
 * The FOM is computed exactly with the BB, using BKZ for the pre-reduction,
 * because the FOM was created with the `ReducerBB` object `red`.
 */
template<typename Int, typename Real>
static double testOneCase(Int m, Int b, IntVec aa, int64_t mindim) {
   int64_t k = aa.length() - 1;
   Int a0inv = computeInva0(b, aa[0]);
   // if (((aa[0] * a0inv) % b) != 1) std::cout << "Error: This inverse is wrong! \n";
   // std::cout << std::dec << "\n";  // We want decimal, not binary.
   if (verbose > 1) {
      Int sumaj = Int(0);  // This will be the sum of coefficients a_j for j > 0.
      for (int64_t j = 1; j <= k; j++)
         sumaj += aa[j];
      std::cout << "Coefficients aa = " << aa << "\n";
      std::cout << " in hexadecimal = ";
      std::cout << " " << "- 0x" << std::hex << conv<uint64_t>(-aa[0]);
      for (int64_t j = 1; j <= k; j++)
         std::cout << " " << "0x" << std::hex << conv<uint64_t>(aa[j]);
      cout << "\n tau = (";
      for (int64_t j = 1; j <= k; j++)
         std::cout << " + 0x" << std::hex << conv<uint64_t>(aa[j]) <<
              " * (__uint128_t) x" << j;
      cout << ") + c;\n";
      cout << "\nInverse of " << aa[0] << " mod b is " << a0inv << " = 0x" << std::hex
            << conv<uint64_t>(a0inv) << endl;
      std::cout << "b - a_k = " << b - aa[k] << "\n";
      std::cout << "b - sum a_j = " << b - sumaj << "\n";
      std::cout << "Modulo m = " << m << "\n";
      std::cout << "Log_2(m) = " << Lg(m) << "\n";
      Int sumaj2 = aa[0] * aa[0];   // Sum of squares of coefficients.
      for (int64_t j = 1; j <= k; j++)
         sumaj2 += aa[j] * aa[j];
      std::cout << "Spectral test bounds:\n";
      std::cout << "a_0^2 + ... + a_k^2 = " << conv<double>(sumaj2) << "\n";
      std::cout << "bound in 2 dim: (b^2 + 1)^{-1/2} = " << 1.0 / sqrt(conv<double>(b * b + 1))
            << "\n";
      std::cout << "bound in k+1 dim: (a_0^2 + ... + a_k^2)^{-1/2} = "
            << 1.0 / sqrt(conv<double>(sumaj2)) << "\n";
      std::cout << "Spectral test in dual lattice:\n\n" << std::dec;
   }
   // Use the LCGLattice object `lcg` to perform the spectral test.
   lcg.setModulus(m);
   lcg.seta(b, maxdim);
   normaDual.computeBounds(-log(m), 1);
   // fomdual.computeMeritSucc(lcg, k+1, maxdim);  // Start in k+1 dim.
   //mwc.setb(b);
   //mwc.setaa(aa, maxdim);   // must implement buildProjection  ****************
   double merit = fomdual.computeMeritSucc(lcg, mindim, maxdim);  // Start in k+1 dim.
   // std::cout << " min merit succ = \n " << merit << "\n";
   //mwc.buildDualBasis(k+2);
   // std::cout << " Initial dual basis:\n " << mwc.getDualBasis() << "\n";
   if (nonsucc) {
      if (verbose > 1) std::cout << "\n";
      proj.setModulus(m);
      merit = min(merit, fomdual.computeMeritNonSucc(lcg, proj));
   }
   return merit;
}

/**
 * We perform a search for good parameters for MWC generators of order `k`, modulus `b = 2^e`,
 * and `numaj` positive coefficients `a_j` for j > 0.
 * The nonzero coefficients `a_j` are generated randomly with `e0` random bits for a_0,
 * `ek` random bits for `a_k`, and `ej` random bits for the other positive coefficients.
 * When `e0 = 0`, we set `a_0 = -1`, otherwise, we have `a_0 > 0$.
 * To impose the condition that `a_k = a_{k-1} = ... = a_{k-c+1}`, just take `numaj = -c`.
 * If `nonsucc`, we look at non-successive projections in the FOM.
 * We try `numMult` random vectors of parameters for a MWC with the given constraints.
 * We test all those for which `m` is a safe prime, and finally we print the best FOM.
 * Some variables such as `maxdim`, `lcg`, `t`, etc., are global instead of passed as parameters.
 */
template<typename Int, typename Real>
static void tryMultipliers(int64_t k, int64_t numaj, int64_t e, int64_t e0, int64_t ek, int64_t ej,
      int64_t numMult, int64_t mindim) {
   Int b = NTL::power(Int(2), e);
   Int m;
   assert(numaj <= k);
   IntVec bestaa;   // Coefficients for the best candidate.
   bestaa.SetLength(k + 1);
   IntVec aa;  // Coefficients for current candidate.
   aa.SetLength(k + 1);
   aa[0] = -1;
   for (int64_t j = 1; j < k; j++)
      aa[j] = 0;
   fomdual.setTVector(t, true);
   numSafe = 0;
   numMaxPer = 0;
   bestmerit = 0.0;

   std::cout << "\n=============================================================\n";
   std::cout << "tryMultipliers: \n";
   std::cout << "Order k = " << std::dec << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   if (e0 == 0) std::cout << "We impose a_0 = -1." << "\n";
   std::cout << "Number of nonzero coefficients a_j besides a_0: " << numaj << "\n";
   std::cout << "Random bits for coefficients, e_0 = " << std::dec << e0 << ",  "
         "e_k = " << ek << ", e_j = " << ej << "\n";
   std::cout << "We try " << numMult << " random vectors of coefficients.\n";
   if (verbose > 1) std::cout << "We found the following candidates: \n\n";
   // fomdual.setVerbosity(4);
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   for (int64_t i = 1; i <= numMult; i++) {
      if (e0 > 0) {
         aa[0] = conv<Int>(LatticeTester::RandBits(e0));
         if ((aa[0] % 2) == 0) aa[0] += 1;  // a_0 must be odd.
      }
      // aa[k] = b - conv<Int>(i);   // If we want a_k very close to b, for testing.
      // aa[k] = b - NTL::power(Int(2), 60) + conv<Int>(LatticeTester::RandBits(58));
      aa[k] = conv<Int>(LatticeTester::RandBits(ek));  // `ek` random bits for a_k.
      for (int64_t j = 1; j < min(numaj, k); j++)
         aa[k-j] = conv<Int>(LatticeTester::RandBits(ej));  // `ej` random bits for a_{k-j}.
      for (int64_t j = 1; j < -numaj; j++)  // To impose equal coefficients, for testing.
         aa[k-j] = aa[k];
      m = computeLCGModulusMWC(b, aa);
      if (m < Int(2)) std::cout << "tryMultipliers: Modulo m < 2 (too small), m = " << m << "\n";
      //if (testOrder (m, b, aa)) {
      bool maxper = mIsSafePrime(m, 100);
      if (maxper) numSafe++;
      else if (!onlySafe) maxper = maxPeriodHalfMWC (m, e);
      if (maxper) {
         numMaxPer++;
         if (verbose > 1) std::cout << "i = " << std::dec << i << ", has max period of p \n";
         merit = testOneCase<Int, Real>(m, b, aa, mindim);
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
   if (verbose > 1) std::cout << "\n\n";
   std::cout << "----------------------------------------------------------------\n";
   std::cout << "Number of safe primes:  " << numSafe << "\n";
   std::cout << "Number with max period: " << numMaxPer << "\n";
   std::cout << "Best FOM: " << bestmerit << "\n";
   std::cout << "Best aa: " << bestaa << "\n\n";
   m = computeLCGModulusMWC(b, bestaa);
   verbose = 2;
   fomdual.setVerbosity(4, 0.01);  // This will print lots of details on the FOM for the winner.
   merit = testOneCase<Int, Real>(m, b, bestaa, mindim);
   std::cout << "Merit = " << merit << "\n";
   /*
   // Now we look at the projections onto `k` coordinates, in the first 12 coordinates.
   NTL::Vec<int64_t> tk; // The t-vector for the FOM.
   tk.SetLength(3);
   tk[0] = tk[1] = tk[2] = 0;
   tk[k-1] = 12;   // Projections onto k coordinates.
   fomdual.setTVector(tk, true);
   merit = min(merit, fomdual.computeMeritNonSucc(lcg, proj));
   std::cout << "Merit for projections onto k coordinates = " << merit << "\n";
   std::cout << "Worst projection onto k coordinates = " << fomdual.getMinMeritProj() << "\n";
   std::cout << "=========================================================================\n\n";
   */

   tmp = clock() - tmp;
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
   verbose = 0;
   fomdual.setVerbosity(2, 0.01);
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: ";
   std::cout << "We make a search for MWC generators of order k, modulo b, \n";
   std::cout << "period (m-1)/2, and apply spectral test (in dual) to retained ones.\n\n";
   fomdual.setVerbosity(2, 0.01);   // Do not print the details for all safe primes.
   verbose = 0;
   onlySafe = true; // If `false`, we look for all MWCs with period p, not only safe primes.
   nonsucc = true;  // We consider projections over non-successive coordinates, specified by t.

   t.SetLength(5);
   t[0] = 12;   // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 12;    // pairs
   t[2] = 10;    // triples
   t[3] = 10;    // etc.
   t[4] = 8;

   int64_t k = 1;   // Order of the MWC.
   // int64_t numaj = 1;
   // int64_t e = 64;    // b = 2^{64}
   // int64_t e0 = 60;   // Number of random bits for a_0. If 0, then a_0 = -1.
   // int64_t ek = 62;   // Number of random bits for a_k.
   // int64_t ej = 58;   // Number of random bits for other a_j's.
   // int64_t maxdim = 12;  // Max dimension for lattice.
   int64_t numMult = 1000 * 1000 * 100;  // 100 millions
   //int64_t numMult = 1000 * 1000;  // one million
   //int64_t numMult = 1000 * 10;  // 10 thousands

   // tryMultipliers<Int, Real>(k, numaj, e, e0, ek, ej, numMult, mindim);

   k = 1;
   tryMultipliers<Int, Real>(k, 1, 64, 0, 62, 0, numMult, k+1);

   k = 2;
   t[1] = 0;  // No pairs when k = 2
   // tryMultipliers<Int, Real>(2, 1, 10, 0, 8, 8, 10, k+1);   // MWC1k2-0-10.res
   // return 0;

   tryMultipliers<Int, Real>(2, 1, 64, 0, 62, 60, numMult, k+1);   // MWC1k2-0-62.res
   tryMultipliers<Int, Real>(2, 2, 64, 0, 62, 60, numMult, k+1);   // MWC1k2-62-60.res
   tryMultipliers<Int, Real>(2, -2, 64, 0, 62, 60, numMult, k+1);   // MWC1k2-62-60.res

   tryMultipliers<Int, Real>(2, 1, 64, 60, 62, 60, numMult, k+1);  // MWCa0k2-0-62.res
   tryMultipliers<Int, Real>(2, 2, 64, 60, 62, 60, numMult, k+1);  // MWCa0k2-62-60.res
   // return 0;

   k = 3;
   t[1] = 0;  // No pairs
   t[2] = 0;  // No triples when k = 3
   tryMultipliers<Int, Real>(3, 1, 64, 0, 62, 60, numMult, k+1);   // MWC1k3-00-62.res
   tryMultipliers<Int, Real>(3, 2, 64, 0, 61, 60, numMult, k+1);   // MWC1k3-0-60-58.res
   tryMultipliers<Int, Real>(3, 3, 64, 0, 61, 58, numMult, k+1);   // MWC1k3-60-58.res
   tryMultipliers<Int, Real>(3, -2, 64, 0, 61, 58, numMult, k+1);   // MWC1k3-60-58.res
   tryMultipliers<Int, Real>(3, -3, 64, 0, 61, 58, numMult, k+1);   // MWC1k3-60-58.res

   tryMultipliers<Int, Real>(3, 1, 64, 60, 62, 58, numMult, k+1);  // MWCa0k3-00-62.res
   tryMultipliers<Int, Real>(3, 2, 64, 60, 60, 58, numMult, k+1);  // MWCa0k3-0-60-58.res
   tryMultipliers<Int, Real>(3, 3, 64, 60, 60, 58, numMult, k+1);  // MWCa0k3-60-58.res

   return 0;
}

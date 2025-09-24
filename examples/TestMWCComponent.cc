/**
 * TestMWCComponent: Testing MWC generators.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/Random.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"

using namespace LatticeTester;
using namespace LatMRG;

template<typename Int>
static void tryMultipliers(int64_t k, int64_t e, int64_t numMult) {
   Int b = getPow2<Int>(e);
   std::cout << "Modulo b = " << b << "\n";
   Int m;
   IntVec aa;
   aa.SetLength(k + 1);
   aa[0] = Int(1);
   for (int64_t i = 0; i < numMult; i++) {
      for (int64_t j = 1; j <= k; j++)
         aa[j] = conv<Int>(LatticeTester::RandBits(30));
      m = computeLCGModulusMWC(b, aa);
      if (mIsPrime(m) <= 1) {
         // std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
         // std::cout << "  vector aa = " << aa << "\n";
         if (mIsSafePrime(m, 50) <= 1) {
            std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
            std::cout << "This modulus is a safe prime!\n";
            IntFactorization<Int> fact(m - Int(1));
            LatMRG::PrimeType ptype = fact.decompToFactorsInv(DECOMP, NULL);
            // std::cout << "  decomp done, num factors: " <<  (m_fact.getFactorList()).size() << "\n";
            if (ptype > 1) {
               std::cout << "  Factorization failed. " << "\n";
            }
            if (isPrimitiveElement(Int(2), fact, m - Int(1)))
               std::cout << "  Also gives a maximal period since 2 is a primitive element.\n";
         }
      }
   }
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: examples of MWC generators.\n\n";

   int64_t k = 3;
   int64_t e = 64;  // b = 2^{64}
   tryMultipliers<Int>(k, e, 100000);
   return 0;

   // Note how large ZZ integers can be initialized in NTL.
   Int m = to_ZZ("67120242110049185336291691984133816319");
   Int m2 = (m - Int(1)) / Int(2);
   // LatMRG::PrimeType ptype = IntFactor<Int>::isPrime(m2, 50);
   std::cout << "\n\n m = " << m << "\n";
   std::cout << "Prime Type of (m-1)/2 : " << IntFactor<Int>::isPrime(m2, 50) << "\n";

   return 0;
}

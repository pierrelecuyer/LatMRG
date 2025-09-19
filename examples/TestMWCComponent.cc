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
static void tryMultipliers(MWCComponent<Int> mwc, int64_t numMult) {
   int64_t k = mwc.getOrder();
   Int m;
   IntVec aa;
   aa.SetLength(k+1);
   aa[0] = Int(1);
   for (int64_t i = 0; i < numMult; i++) {
      for (int64_t j = 1; j <= k; j++)
         aa[j] = conv<Int>(LatticeTester::RandBits(30));
      mwc.setaa(aa);
      m = mwc.computeLCGModulus();
      // std::cout << "Modulus m = " << mwc.getLCGm() << "\n";
      if (mwc.mIsPrime()) {
         std::cout << "i = " << i << "\n";
         std::cout << "Modulo m = " << m << "  is prime.\n";
         std::cout << "  vector aa = " << aa << "\n";
         if (mwc.maxPeriod1Pow2(100)) {
            std::cout << "  Also gives a maximal period.\n";
            if (mwc.mIsSafePrimePow2(100))
               std::cout << "  m is also a safe prime.\n";
         }
      }
   }
}



int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: examples of MWC generators.\n\n";

   int64_t k = 3;
   MWCComponent<Int> mwc(k);
   mwc.setbPow2(32);
   std::cout << "Modulo b = " << mwc.getb() << "\n";
   // tryMultipliers(mwc, 10000);

   IntVec aa;
   aa.SetLength(k+1);
   aa[0] = Int(1);
   aa[1] = Int(24585293);
   aa[2] = Int(206700418);
   aa[3] = Int(4768041);
   mwc.setaa(aa);
   Int m = mwc.computeLCGModulus();
   std::cout << "Modulo m = " << m << "\n";
   if (mwc.maxPeriod1Pow2(100)) {
      std::cout << "  Gives a maximal period.\n";
      if (mwc.mIsSafePrimePow2(100))
         std::cout << "  m is also a safe prime.\n";
   }

   return 0;
}

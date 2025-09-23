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
   aa.SetLength(k + 1);
   aa[0] = Int(1);
   for (int64_t i = 0; i < numMult; i++) {
      for (int64_t j = 1; j <= k; j++)
         aa[j] = conv<Int>(LatticeTester::RandBits(30));
      mwc.setaa(aa);
      m = mwc.computeLCGModulus();
      if (mwc.mIsPrime() <= 1) {
         // std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
         // std::cout << "  vector aa = " << aa << "\n";
         if (mwc.mIsSafePrimePow2(50) <= 1) {
         //if (mwc.maxPeriod1Pow2(100)) {
            //if (mwc.mIsSafePrimePow2(50) > 1) {
               std::cout << "i = " << i << ", Modulo m = " << m << "  is prime.\n";
               std::cout << "This modulus is a safe prime!\n";
                // std::cout << "  Also gives a maximal period.\n";
                // if (mwc.mIsSafePrimePow2(100))
                // std::cout << "  m is a safe prime.\n";
                IntFactorization<Int> m_fact(m - Int(1));
                // std::cout << "  m_fact created " << "\n";
                LatMRG::PrimeType ptype = m_fact.decompToFactorsInv(DECOMP, NULL);
                // std::cout << "  decomp done, num factors: " <<  (m_fact.getFactorList()).size() << "\n";
                if (ptype > 1) {
                std::cout << "  Factorization failed. " << "\n";
                // return;
                }
                if (isPrimitiveElement(Int(2), m_fact, m - Int(1)))
                std::cout << "  Also gives a maximal period since 2 is a primitive element.\n";
//            }
         }
      }
   }
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCComponent: examples of MWC generators.\n\n";

   int64_t k = 3;
   MWCComponent<Int> mwc(k);
   mwc.setbPow2(64);
   std::cout << "Modulo b = " << mwc.getb() << "\n";
   tryMultipliers(mwc, 100000);
   return 0;

   // Note how large ZZ integers can be initialized in NTL.
   Int m = to_ZZ("67120242110049185336291691984133816319");
   Int m2 = (m - Int(1)) / Int(2);
   // LatMRG::PrimeType ptype = IntFactor<Int>::isPrime(m2, 50);
   std::cout << "\n\n m = " << m << "\n";
   std::cout << "Prime Type of (m-1)/2 : " << IntFactor<Int>::isPrime(m2, 50) << "\n";

   /*
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
    */
   return 0;
}

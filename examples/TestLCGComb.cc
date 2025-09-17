/**
 * TestLCGComb: Examples of combined LCGs.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/MRGCombined.h"
#include "latmrg/LCGCombined.h"

using namespace LatMRG;

template<typename Int>
static void computeCombinedLCG(int64_t numcomp, Int mm[], Int aa[]) {
   LCGCombined<Int> lcgcomb(numcomp, mm, aa);
   lcgcomb.computeCombination();
   std::cout << "Parameters of the combined LCG:\n";
   std::cout << "Modulo m = " << lcgcomb.getModulus() << "\n";
   std::cout << "Multiplier a = " << lcgcomb.geta() << "\n\n";
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestCombLCG: examples of combined LCGs.\n\n";
   //LCGComponent<Int> lcg1(m), lcg2(m), lcg3(m), lcg4(m);  // The components.
   //LCGCombined<Int> lcgcomb;   // The combination.

   // Combined 32-bit LCG from L'Ecuyer (1986, 1988).
   std::cout << "Two-component combined 32-bit LCG from L'Ecuyer (1986, 1988).\n";
   Int m86[2] = { Int(2147483563), Int(2147483399) };
   Int a86[2] = { Int(40014), Int(40692) };
   computeCombinedLCG(2, m86, a86);

   // Combined 32-bit LCG from L'Ecuyer and Tezuka (1991).
   std::cout << "Two-component combined 32-bit LCG from L'Ecuyer and Tezuka (1991).\n";
   Int m91[2] = { Int(2147483647), Int(2145483479) };
   Int a91[2] = { Int(26756), Int(30318) };
   computeCombinedLCG(2, m86, a86);

   // Combined LCG from L'Ecuyer and Andres (1997).
   std::cout << "Four-component combined 32-bit LCG from L'Ecuyer and Andres (1997).\n";
   Int m97[4] = { Int(2147483647), Int(2147483543), Int(2147483423), Int(2147483323) };
   Int a97[4] = { Int(45991), Int(207707), Int(138556), Int(49689) };
   computeCombinedLCG(4, m86, a86);

   return 0;
}

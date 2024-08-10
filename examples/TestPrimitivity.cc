/**
 * This example tests the functions that deal with integer factorization, primality testing,
 * finding certain types of prime numbers, and testing for the primitivity of an integer
 * modulo m or of a polynomial modulo another polynomial f(x).
 **/

// The code to define the Int and Real types.  Here we must recompile to change it.
//#define TYPES_CODE  LD     // Int = int64_t, Real = double
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
// #include "latticetester/ReducerBB.h"
#include "latmrg/FlexModInt.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"
// #include "latmrg/PrimesFinder.h"

#include "latmrg/LCGComponent.h"
#include "latmrg/MRGComponent.h"


using namespace LatMRG;
// using namespace LatticeTester;

template<typename Int>
static void TestPrimitiveInt(const Int &m, const Int &a) {
   std::cout << "TestPrimitiveInt for m = " << m << ", a = " << a << "\n";
   IntFactorization<Int> fac(m-1);
   fac.factorizePlus();
   std::cout << fac.toString();
   bool isP = isPrimitiveElement(a, fac, m);
   std::cout << a << " is a primitive element mod " << m << "? " << isP << "\n\n";
}

template<typename Int>
static void TestPeriodLCG (const Int &m, const Int &a, DecompType decompm1, const char *filem1) {
   LCGComponent<Int> lcg(m, decompm1, filem1);
   std::cout << "TestPeriodLCG for m = " << m << ", a = " << a << "\n";
   if (lcg.maxPeriod(a)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}

template<typename Int>
static void TestPeriodLCG (const Int &b, int e, const Int &c, const Int &a, DecompType decompm1, const char *filem1) {
   LCGComponent<Int> lcg(b, e, c, decompm1, filem1);
   std::cout << "TestPeriodLCG for m = " << lcg.getModulus() << ", a = " << a << "\n";
   if (lcg.maxPeriod(a)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}

template<typename Int>
static void TestPeriodMRG (const Int &m, int k, const NTL::vector<Int> &aa, DecompType decompm1, const char *filem1, DecompType decompr,
      const char *filer) {
   std::cout << "TestPeriodMRG for m = " << m << ", aa = " << aa << "\n";
   MRGComponent<Int> mrg(m, k, decompm1, filem1, decompr, filer);
   std::cout << "TestPeriodMRG for m = " << mrg.getModulus() << ", aa = " << aa << "\n";
   if (mrg.maxPeriod(aa)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}



int main() {
   typedef NTL::ZZ Int;
   std::cout << "/nTestPrimitivity with Types: " << strFlexTypes << "\n\n";

   Int m, a, r, b, c;
   int e;
   DecompType decompm1 = DECOMP;
   char *filem1 = NULL;
   // DecompType decompr = DECOMP;
   // char *filer = NULL;

   m = 2147483647;
   a = 16807;        // An LCG multiplier
   b = 2;  e = 31;  c = -1;

   TestPrimitiveInt(m, a);
   TestPeriodLCG (m, a, decompm1, filem1);
   TestPeriodLCG (b, e, c, a, DECOMP, NULL);

   int k = 3;
   Int mm(9223372036854773561);
   NTL::vector<Int> aa;
   aa.SetLength(k+1);
   aa[0] = 1;
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   // std::cout << "Before TestPeriodMRG, aa = " << aa << "\n";
   TestPeriodMRG (mm, 3, aa, DECOMP, NULL, DECOMP_PRIME, NULL);

   /*
   mm = 101;
   setModulusIntP<Int>(mm);
   r = (mm * mm * mm - 1) / (mm - 1);
   aa[0] = 1;
   aa[1] = 28;
   aa[2] = 49;
   aa[3] = 8;
   std::cout << "\nTesting for a primitive polynomial with m = " << mm << ", aa = " << aa << ", r = " << r << "\n";
   IntFactorization<Int> fm (mm-1);
   fm.factorizePlus();
   IntFactorization<Int> fr (r);
   fr.setStatus(PRIME);
   isP = isPrimitiveElement(aa[3], fm, mm);
   std::cout << aa[3] << " is a primitive element mod " << mm << "? " << isP << "\n";
   isP = isPrimitive<Int>(aa, mm, fm, fr);
   std::cout << "Is it a primitive polynomial?  " << isP << "\n";

   mm = 9223372036854773561;
   setModulusIntP<Int>(mm);
   r = (mm * mm * mm - 1) / (mm - 1);
   aa[0] = 1;
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   std::cout << "\nTesting for a primitive polynomial with m = " << mm << ", aa = " << aa << ", r = " << r << "\n";
   fm.setNumber (mm-1);
   fm.factorizePlus();
   fr.setNumber (r);
   fr.setStatus(PRIME);
   isP = isPrimitiveElement(aa[3], fm, mm);
   std::cout << aa[3] << " is a primitive element mod " << mm << "? " << isP << "\n";
   isP = isPrimitive<Int>(aa, mm, fm, fr);
   std::cout << "Is it a primitive polynomial?  " << isP << "\n";
*/

   return 0;
}


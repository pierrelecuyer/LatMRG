/**
 * This example tests the functions that deal with integer factorization, primality testing,
 * finding certain types of prime numbers, and testing for the primitivity of an integer
 * modulo m or of a polynomial modulo another polynomial f(x).
 * For his, it computes the period length of certain LCGs and MRGs.
 **/

// The code to define the Int and Real types.  Here we must recompile to change it.
//#define TYPES_CODE  LD     // Int = int64_t, Real = double
// #define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latmrg/FlexModInt.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"
// #include "latmrg/PrimesFinder.h"

#include "latmrg/LCGComponent.h"
#include "latmrg/MRGComponent.h"

using namespace LatMRG;

template<typename Int>
static void TestPrimitiveInt(const Int &m, const Int &a) {
   std::cout << "TestPrimitiveInt for m = " << m << ", a = " << a << "\n";
   IntFactorization<Int> fac(m-1);
   fac.factorizePlus();
   std::cout << "Factorization of " << fac.toString();
   bool isP = isPrimitiveElement(a, fac, m);
   std::cout << a << " is a primitive element mod " << m << "? " << isP << "\n\n";
}

template<typename Int>
static void TestPeriodLCG (const Int &m, const Int &a, DecompType decomp, const char *filename, bool increment) {
   LCGComponent<Int> lcg(m, decomp, filename, increment);
   std::cout << "TestPeriodLCG for m = " << m << ", a = " << a << ", increment = " << increment << "\n";
   if (lcg.maxPeriod(a)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}

template<typename Int>
static void TestPeriodLCG (const Int &b, int e, const Int &c, const Int &a, DecompType decomp, const char *file, bool increment) {
   LCGComponent<Int> lcg(b, e, c, decomp, file, increment);
   std::cout << "TestPeriodLCG for m = " << lcg.getModulus() << ", a = " << a << ", increment = " << increment << "\n";
   if (lcg.maxPeriod(a)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}

template<typename Int>
static void TestPeriodMRG (const Int &m, int k, const NTL::Vec<Int> &aa, DecompType decompm1, const char *filem1, DecompType decompr,
      const char *filer) {
   MRGComponent<Int> mrg(m, k, decompm1, filem1, decompr, filer);
   std::cout << "TestPeriodMRG for m = " << mrg.getModulus() << ", aa = " << aa << "\n";
   if (mrg.maxPeriod(aa)) std::cout << "Max period: YES \n\n";
   else std::cout << "Max period: NO \n\n";
}

int main() {
   typedef NTL::ZZ Int;
   std::cout << "\nTestPeriod with Int = ZZ  \n\n";
   Int m, a, r, b, c;
   int e;
   char filem1[] = "decomp2147483646.txt";

   m = 2147483647;
   a = 16807;        // An LCG multiplier
   b = 2;  e = 31;  c = -1;

   TestPrimitiveInt<Int>(m, a);
   TestPeriodLCG<Int> (m, a, DECOMP_WRITE, filem1, false);
   std::cout << "Creating an IntFactorization from filem1 = " << filem1 << "\n";
   IntFactorization<Int> fact = IntFactorization<Int>(filem1);
   std::cout << fact.toString();
   TestPeriodLCG<Int> (b, e, c, a, DECOMP_READ, filem1, false);

   TestPeriodLCG<Int> (m, a, DECOMP, NULL, true);
   TestPeriodLCG<Int> (m+1, (Int)16801, DECOMP, NULL, true);

   int k = 3;
   NTL::Vec<Int> aa;
   aa.SetLength(k+1);
   // This one works only with Int = ZZ, numbers are too large for int64_t.
   Int mm(9223372036854773561);
   aa[0] = 1;
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   TestPeriodMRG<Int> (mm, 3, aa, DECOMP, NULL, DECOMP_PRIME, NULL);

   return 0;
}


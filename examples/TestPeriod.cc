/**
 * This example tests the functions that deal with integer factorization, testing for the
 * primitivity of an integer modulo m or of a polynomial modulo another polynomial f(x),
 * and testing if an LCG or an MRG has maximal period.
 */

// The code to define the Int and Real types.  Here we must recompile to change it.
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MRGComponent.h"

using namespace LatMRG;

std::string boolYN[2] = { "NO", "YES" };

// Assumes that `m` is prime and tests if `a` is primitive mod m.
// Performs the required factorization of `m-1` internally.
template<typename Int>
static void TestPrimitiveInt(const Int &m, const Int &a) {
   std::cout << "TestPrimitiveInt for m = " << m << ", a = " << a << "\n";
   IntFactorization<Int> fac(m-1);
   fac.factorizePlus();
   std::cout << "Factorization of " << fac.toString();
   bool isP = isPrimitiveElement(a, fac, m);
   std::cout << a << " is a primitive element mod " << m << "? " << isP << "\n\n";
}

// Tests if this LCG has full period by calling `maxPeriod` component.
// When `m` is prime, this is equivalent to test if `a` is primitive mod m.
// The required factorization can be given in a file.
template<typename Int>
static void TestPeriodLCG (const Int &m, const Int &a, DecompType decomp, const char *filename, bool increment) {
   LCGComponent<Int> lcg(m, decomp, filename, increment);
   std::cout << "TestPeriodLCG for m = " << m << ", a = " << a << ", increment = " << boolYN[increment] << "\n";
   std::cout << "Max period: " << boolYN[lcg.maxPeriod(a)] << "\n\n";
}

// Same as previous function, but `m` is given in a different format.
template<typename Int>
static void TestPeriodLCG (const Int &b, int e, const Int &c, const Int &a, DecompType decomp, const char *file, bool increment) {
   LCGComponent<Int> lcg(b, e, c, decomp, file, increment);
   std::cout << "TestPeriodLCG for m = " << lcg.getModulus() << ", a = " << a << ", increment = " << boolYN[increment] << "\n";
   std::cout << "Max period: " << boolYN[lcg.maxPeriod(a)] << "\n\n";
}

// Constructs an `MRGComponent` object and tests if it has maximal period.
template<typename Int>
static void TestPeriodMRG (const Int &m, int k, const NTL::Vec<Int> &aa, DecompType decompm1, const char *filem1, DecompType decompr,
      const char *filer) {
   MRGComponent<Int> mrg(m, k, decompm1, filem1, decompr, filer);
   std::cout << "TestPeriodMRG for m = " << mrg.getModulus() << ", aa = " << aa << "\n";
   bool maxper = mrg.maxPeriod(aa);
   std::cout << "Max period: " << boolYN[maxper] << "\n\n";
}

int main() {
   typedef NTL::ZZ Int;
   std::cout << "\nTestPeriod with Int = ZZ  \n\n";
   Int m, a, r, b, c;
   int64_t e;
   char filem1[] = "decomp2147483646.txt";

   m = 2147483647;
   a = 16807;        // An LCG multiplier
   b = 2;  e = 31;  c = -1;

   // We factorize m-1, print the factors, and put them in a file.
   std::cout << "Creating an IntFactorization for m-1 and write the factors in filem1 = " << filem1 << "\n";
   IntFactorization<Int> factm1(m-1);
   factm1.factorize();
   std::cout << factm1.toString();
   std::ofstream fout(filem1);
   fout << factm1.toString();

   // Here we test if a is a primitive element mod m. The factorizations are done internally.
   TestPrimitiveInt<Int>(m, a);

   // We test if this LCG has full period and save the decomposition of m-1 in a file.
   TestPeriodLCG<Int> (m, a, DECOMP_WRITE, filem1, false);
   std::cout << "Creating an IntFactorization from filem1 = " << filem1 << "\n";
   IntFactorization<Int> fact(filem1);
   std::cout << fact.toString();

   // We test the same thing again, but this time we read the decomposition from the file.
   TestPeriodLCG<Int> (b, e, c, a, DECOMP_READ, filem1, false);

   // Same test again, but the decomposition is now done internally and not saved.
   TestPeriodLCG<Int> (m, a, DECOMP, NULL, true);

   // Here we test the period for the LCG with modulus 2^{31}.
   TestPeriodLCG<Int> (m+1, (Int)16801, DECOMP, NULL, true);

   // We now test the period for an MRG of order k=3 and larger modulus.
   // The required factorizations are performed internally.
   // This one works only with Int = ZZ, the numbers are too large for int64_t.
   int k = 3;
   NTL::Vec<Int> aa;
   aa.SetLength(k+1);
   Int mm(9223372036854773561);
   aa[0] = 1;
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   TestPeriodMRG<Int> (mm, 3, aa, DECOMP, NULL, DECOMP_PRIME, NULL);

   return 0;
}


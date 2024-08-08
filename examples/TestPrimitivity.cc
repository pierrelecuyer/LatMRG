/**
 * This example tests the functions that deal with integer factorization, primality testing,
 * finding certain types of prime numbers, and testing for the primitivity of an integer
 * modulo m or of a polynomial modulo another polynomial f(x).
 **/

// The code to define the Int and Real types.  Here we must recompile to change it.
#define TYPES_CODE  LD     // Int = int64_t, Real = double
//#define TYPES_CODE  ZD     // Int = ZZ, Real = double

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
#include "latmrg/PrimitiveInt.h"
#include "latmrg/PrimitivePoly.h"
// #include "latmrg/PrimesFinder.h"

using namespace LatMRG;

Int m(2147483647);   // Prime modulus m
Int a(16807);        // An LCG multiplier
bool isP;

NTL::ZZ mm, r;
NTL::vector<NTL::ZZ> aa(4); // Vector a has size 4, a[j] contains a_j.
// typename ModInt<NTL::ZZ>::PolX f;

int main() {
   std::cout << "Types: " << strFlexTypes << "\n";
   std::cout << "TestPrimitivity modulo m = " << m << "\n\n";
   IntFactorization<Int> fac(m - 1);
   fac.factorize();
   std::cout << fac.toString();
   fac.calcInvFactors();
   isP = isPrimitiveElement(a, fac, m);
   std::cout << a << " is a primitive element mod " << m << "? " << isP << "\n";
   a = 123456;
   isP = isPrimitiveElement(a, fac, m);
   std::cout << a << " is a primitive element mod " << m << "? " << isP << "\n";

//   mm = 9223372036854773561;
//   aa[0] = 1;
//   aa[1] = 1145902849652723;
//   aa[2] = 0;
//   aa[3] = -1184153554609676;

   mm = 101;
   setModulus<NTL::ZZ>(mm);
   // ModInt<Int>::IntP::init(mm);        // Sets the modulus to m.
   r = (mm * mm * mm - 1) / (mm - 1);
   aa[0] = 1;
   aa[1] = 28;
   aa[2] = 49;
   aa[3] = 8;
   std::cout << "\nTesting for a primitive polynomial with m = " << mm << ", aa = " << aa << ", r = " << r << "\n";
   IntFactorization<NTL::ZZ> fm (mm-1);
   fm.factorizePlus();
   IntFactorization<NTL::ZZ> fr (r);
   fr.setStatus(PRIME);
   isP = isPrimitiveElement(aa[3], fm, mm);
   std::cout << aa[3] << " is a primitive element mod " << mm << "? " << isP << "\n";

   isP = isPrimitive<NTL::ZZ>(aa, mm, fm, fr);
   std::cout << "Is it a primitive polynomial?  " << isP << "\n";

   mm = 9223372036854773561;
   setModulus<NTL::ZZ>(mm);
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

   isP = isPrimitive<NTL::ZZ>(aa, mm, fm, fr);
   std::cout << "Is it a primitive polynomial?  " << isP << "\n";

   return 0;
}


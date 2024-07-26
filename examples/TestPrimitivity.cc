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
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/PrimitiveInt.h"
#include "latmrg/PrimitivePoly.h"
// #include "latmrg/PrimesFinder.h"

using namespace LatMRG;

Int m(101);      // Modulus m = 101
Int a(34);       // An LCG multiplier
bool isP;



int main() {
   std::cout << "Types: " << strFlexTypes << "\n";
   std::cout << "TestPrimitivity \n\n";
   IntFactorization<Int> f(m-1);
   f.factorize();
   f.calcInvFactors();
   isP = isPrimitiveElement(a, f, m);
   std::cout
         << a << " is a primitive element mod " << m <<"? " << isP << "\n\n";
   return 0;
}


/**
 * This example tests the functions that deal with integer factorization, primality testing,
 * finding certain types of prime numbers, and testing for the primitivity of an integer
 * modulo m or the primitivity of a polynomial f(x) modulo m.
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
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/PrimitiveInt.h"
#include "latmrg/PrimitivePoly.h"
// #include "latmrg/PrimesFinder.h"

using namespace LatMRG;

Int m(101);      // Modulus m = 101
Int a(33);       // An LCG multiplier



int main() {
   std::cout << "Types: " << strFlexTypes << "\n";
   std::cout << "TestPrimitivity \n\n";

   std::cout
         << "We see that the dual of the projection differs from the projection of the dual! \n\n";
   return 0;
}


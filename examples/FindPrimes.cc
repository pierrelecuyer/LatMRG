// This defines the Int type. We must recompile to change it.
#define TYPES_CODE  ZD     // Int = ZZ

#include "latmrg/PrimesFinder.h"

using namespace LatticeTester;
using namespace LatMRG;

int main() {
   // Just find prime numbers in some range.
   findPrime<Int>(31, 2, false, cout);
   findPrime<Int>(61, 2, false, cout);
   findPrime<Int>(1, 61, 2, false, false, cout);
   findPrime<Int>(1, 61, -100, -1, false, false, cout);
   findPrime<Int>(63, 4, true, cout);

   // Find special prime numbers m, for which `(m-1)/2` and `(m^k-1)/(m-1)` are also prime.
   findPrime<Int>(5, 63, 3, true, false, cout);
   cout << "The three values should be the same as in Table 1 of `rLEC99b`.\n";
   return 0;
}

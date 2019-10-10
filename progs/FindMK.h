#ifndef LATMRG_FINDMK_H
#define LATMRG_FINDMK_H
#include "latmrg/Primes.h"
extern std::ostream* out;

namespace LatMRG {
  template<typename Int> struct MKSearch {
    int k, e, c1 = 0, c2;
    bool power = false, safe, facto = false;

    int FindMK() {
      // trouver 3 modules m proches de 2^31 pour des MRGs d'ordre k
      Primes<Int> primes;
      if (power) primes.find (k, e, c1, safe, facto, *out);
      else primes.find (k, e, c1, c2, safe, facto, *out);
      return 0;
    }
  };
}

#endif

#include "latmrg/Primes.h"

template<typename Int> struct MKSearch {
  int k, e, c1 = 0, c2;
  bool power = false, safe, facto;

  int FindMK() {
    std::ofstream fout ("default.res");
    // trouver 3 modules m proches de 2^31 pour des MRGs d'ordre k
    Primes<Int> primes;
    if (power) primes.find (k, e, c1, safe, facto, fout);
    else primes.find (k, e, c1, c2, safe, facto, fout);
    return 0;
  }
};

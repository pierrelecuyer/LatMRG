#include "latmrg/LatTestPalpha.h"

namespace LatMRG {
  template class LatTestPalpha<std::int64_t, double>;
  template class LatTestPalpha<NTL::ZZ, double>;
  template class LatTestPalpha<NTL::ZZ, NTL::RR>;
}

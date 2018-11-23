#include "latmrg/LatTestSpectral.h"

namespace LatMRG {
  template class LatTestSpectral<std::int64_t, double>;
  template class LatTestSpectral<NTL::ZZ, double>;
  template class LatTestSpectral<NTL::ZZ, NTL::RR>;
}

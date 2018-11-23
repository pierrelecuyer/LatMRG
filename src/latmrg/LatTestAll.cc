#include "latmrg/LatTestAll.h"

namespace LatMRG {
  template class LatTestAll<std::int64_t, double>;
  template class LatTestAll<NTL::ZZ, double>;
  template class LatTestAll<NTL::ZZ, NTL::RR>;
}

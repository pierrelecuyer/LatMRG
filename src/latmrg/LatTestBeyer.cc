#include "latmrg/LatTestBeyer.h"

namespace LatMRG {
  template class LatTestBeyer<std::int64_t, double>;
  template class LatTestBeyer<NTL::ZZ, double>;
  template class LatTestBeyer<NTL::ZZ, NTL::RR>;
}

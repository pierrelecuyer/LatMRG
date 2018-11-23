#include "latmrg/MRGLatticeLac.h"

namespace LatMRG {
  template class MRGLatticeLac<std::int64_t, double>;
  template class MRGLatticeLac<NTL::ZZ, double>;
  template class MRGLatticeLac<NTL::ZZ, NTL::RR>;
}

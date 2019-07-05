#include "latmrg/LCGCarryLattice.h"

namespace LatMRG {
  template class LCGCarryLattice<std::int64_t, double>;
  template class LCGCarryLattice<NTL::ZZ, double>;
  template class LCGCarryLattice<NTL::ZZ, NTL::RR>;
}

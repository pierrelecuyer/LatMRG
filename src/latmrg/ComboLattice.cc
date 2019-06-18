#include "latmrg/ComboLattice.h"

namespace LatMRG {
  template class ComboLattice<std::int64_t, double>;
  template class ComboLattice<NTL::ZZ, double>;
  template class ComboLattice<NTL::ZZ, NTL::RR>;
}

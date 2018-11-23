#include "latmrg/MMRGLattice.h"

namespace LatMRG {
  template class MMRGLattice<std::int64_t, double>;
  template class MMRGLattice<NTL::ZZ, double>;
  template class MMRGLattice<NTL::ZZ, NTL::RR>;
}

#include "latmrg/MRGLatticeFactory.h"

namespace LatMRG {
  template class MRGLatticeFactory<std::int64_t, double>;
  template class MRGLatticeFactory<NTL::ZZ, double>;
  template class MRGLatticeFactory<NTL::ZZ, NTL::RR>;
}

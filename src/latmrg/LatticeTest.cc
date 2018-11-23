#include "latmrg/LatticeTest.h"

namespace LatMRG {
  template class LatticeTest<std::int64_t, double>;
  template class LatticeTest<NTL::ZZ, double>;
  template class LatticeTest<NTL::ZZ, NTL::RR>;
}

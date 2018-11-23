#include "latmrg/TestProjections.h"

namespace LatMRG{
  template class TestProjections<std::int64_t, double>;
  template class TestProjections<NTL::ZZ, double>;
  template class TestProjections<NTL::ZZ, NTL::RR>;
}

#include "latmrg/mrgtypes/MixmaxMMRG.h"

#include <NTL/ZZ.h>

namespace LatMRG
{
  template class MixmaxMMRG<std::int64_t>;
  template class MixmaxMMRG<NTL::ZZ>;
} // end namespace

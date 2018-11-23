#include "latmrg/ParamReaderExt.h"

namespace LatMRG {
  template class ParamReaderExt<std::int64_t, double>;
  template class ParamReaderExt<NTL::ZZ, double>;
  template class ParamReaderExt<NTL::ZZ, NTL::RR>;
}

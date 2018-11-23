#include "latmrg/ReportLat.h"

namespace LatMRG {
  template class ReportLat<std::int64_t, double>;
  template class ReportLat<NTL::ZZ, double>;
  template class ReportLat<NTL::ZZ, NTL::RR>;
}

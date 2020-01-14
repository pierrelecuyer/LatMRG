#include "latmrg/Test.h"

namespace LatMRG {
  namespace Reductions {
    template void reduceFull(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat, std::int64_t);
    template void reduceFull(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat, std::int64_t);
    template void reduceFull(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat, std::int64_t);

    template void reduceBKZ(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    template void reduceBKZ(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    template void reduceBKZ(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);

    template void reduceLLL(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    template void reduceLLL(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    template void reduceLLL(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);

    template void reduceMink(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    template void reduceMink(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    template void reduceMink(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);
  }
}

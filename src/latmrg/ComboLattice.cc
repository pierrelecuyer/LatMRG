#include "latmrg/ComboLattice.h"

namespace LatMRG {
  // Instatiations of the temaplates
  template MRGLattice<std::int64_t, double>* getLatCombo(std::vector<MRGComponent<std::int64_t>*>&, int);
  template MRGLattice<NTL::ZZ, double>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>*>&, int);
  template MRGLattice<NTL::ZZ, NTL::RR>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>*>&, int);

  template class ComboLattice<std::int64_t, double>;
  template class ComboLattice<NTL::ZZ, double>;
  template class ComboLattice<NTL::ZZ, NTL::RR>;
} // end namespace LatMRG

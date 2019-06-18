#ifndef LATMRG_COMBOLATTICE_H
#define LATMRG_COMBOLATTICE_H

#include <vector>

#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

  /**
   * This class represents a combined MRG.
   * It stores a vector of `MRGComponent` and computes the equivalent MRG to
   * their combination. Note that MRGComponenets do not have to be MRG
   * themselves. This only overrides the `toString()` method and provides a
   * constructor that populates the MRGLattice super-class to use this class in
   * methods without an overload specific to it.
   * */
  template<typename Int, typename Dbl>
    class ComboLattice: public MRGLattice<Int, Dbl> {
      public:
        typedef Dbl Float;
        typedef Int Integ;
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;

        ComboLattice(std::vector<MRGComponent<Int>>& comp);

        std::string toString() const override;
        
      private:
        std::vector<MRGComponent<Int>> m_comp;
    };

  extern template class ComboLattice<std::int64_t, double>;
  extern template class ComboLattice<NTL::ZZ, double>;
  extern template class ComboLattice<NTL::ZZ, NTL::RR>;
}
#endif

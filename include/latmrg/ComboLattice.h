#ifndef LATMRG_COMBOLATTICE_H
#define LATMRG_COMBOLATTICE_H

#include <vector>

#include "latticetester/Util.h"

#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

  /**
   * This function computes and returns a `MRGLattice` that is the lattice of the
   * combination of generators described in the `MRGComponent`s in the vector `comp`.
   * This dynamically allocates memory to returned pointer. and needs to be deleted
   * afterwards.
   * */
  template<typename Int, typename Dbl>
    MRGLattice<Int, Dbl>* getLatCombo(std::vector<MRGComponent<Int>>& comp, int maxDim) {
      typedef NTL::vector<Int> IntVec;
      int num_comp = comp.size();
      int k = 0;
      Int modulo = Int(1);
      for (int i = 0; i < num_comp; i++) {
        modulo *= comp[i].getM();
        k = std::max(k, comp[i].k);
      }
      // Filling up vector A
      IntVec A(k+1);
      for (int j = 0; j < num_comp; j++) {
        // Assumes the modulus of the components are relatively primes
        Int n;
        {
          Int b, c, d, e;
          LatticeTester::Euclide(modulo/comp[j].getM(), comp[j].getM(), n, b, c, d, e);
        }
        n %= comp[j].getM();

        for (int i = 1; i <= k; i++) {
          A[i] += comp[j].a[i]*n*modulo/comp[j].getM();
        }
      }
      // Modulo only once
      for (int i = 1; i <= k; i++) {
        A[i] %= modulo;
      }

      return new MRGLattice<Int, Dbl>(modulo, A, maxDim, k, FULL);
    }

  extern template MRGLattice<std::int64_t, double>* getLatCombo(std::vector<MRGComponent<std::int64_t>>& comp, int maxDim);
  extern template MRGLattice<NTL::ZZ, double>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>>& comp, int maxDim);
  extern template MRGLattice<NTL::ZZ, NTL::RR>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>>& comp, int maxDim);

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

        /**
         * Creates a ComboLattice for the set of MRG described in `comp` with a
         * MRG lattice as described in `lat`.
         *
         * Ideally, `lat` has been initialized as a MRG
         * */
        ComboLattice(std::vector<MRGComponent<Int>>& comp,
            MRGLattice<Int, Dbl>& lat);

        /**
         * Prints a string describing all the componenets and the equivalent MRG.
         * */
        std::string toString() const override;

      private:
        std::vector<MRGComponent<Int>> m_comp;

        /**
         * The number of components stored in this object.
         * */
        int m_number;
    }; // end class ComboLattice

  template<typename Int, typename Dbl>
    ComboLattice<Int, Dbl>::ComboLattice(std::vector<MRGComponent<Int>>& comp,
        MRGLattice<Int, Dbl>& lat) : MRGLattice<Int, Dbl>(lat){
      m_comp = comp;
      m_number = comp.size();
    }

  extern template class ComboLattice<std::int64_t, double>;
  extern template class ComboLattice<NTL::ZZ, double>;
  extern template class ComboLattice<NTL::ZZ, NTL::RR>;

} // end namespace LatMRG
#endif

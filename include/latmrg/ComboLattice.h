#ifndef LATMRG_COMBOLATTICE_H
#define LATMRG_COMBOLATTICE_H

#include <vector>

#include "latticetester/Util.h"

#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"
#include "latmrg/MWCLattice.h"

namespace LatMRG {

  /**
   * This function computes and returns a `MRGLattice` that is the lattice of the
   * combination of generators described in the `MRGComponent`s in the vector `comp`.
   * This dynamically allocates memory to returned pointer. and needs to be deleted
   * afterwards.
   * */
  template<typename Int, typename Dbl>
    MRGLattice<Int, Dbl>* getLatCombo(std::vector<MRGComponent<Int>*>& comp, int maxDim) {
      typedef NTL::vector<Int> IntVec;
      int num_comp = comp.size();
      int k = 0;
      Int modulo = Int(1);
      for (int i = 0; i < num_comp; i++) {
        Int mod_gen(1);
        if (comp[i]->get_type() == MRG) {
          k = std::max(k, comp[i]->getK());
          mod_gen = comp[i]->getM();
        } else if (comp[i]->get_type() == MWC) {
          mod_gen = MWCEquiv::LCGMod(comp[i]->m_MWCb, comp[i]->getA());
          k = std::max(k, 1);
        }
        modulo *= mod_gen;
      }
      // Filling up vector A
      IntVec A(k+1);
      for (int j = 0; j < num_comp; j++) {
        // Assumes the modulus of the components are relatively primes
        Int mod_gen(1);
        if (comp[j]->get_type() == MRG) {
          mod_gen = comp[j]->getM();
        } else if (comp[j]->get_type() == MWC) {
          mod_gen = MWCEquiv::LCGMod(comp[j]->m_MWCb, comp[j]->getA());
        }
        Int n;
        {
          Int b, c, d, e;
          LatticeTester::Euclide(modulo/mod_gen, mod_gen, n, b, c, d, e);
        }
        n %= mod_gen;

        IntVec a;
        if (comp[j]->get_type() == MRG) {
          a = comp[j]->getA();
        } else if (comp[j]->get_type() == MWC) {
          a.resize(1);
          a[0] = MWCEquiv::LCGCoeff(comp[j]->m_MWCb, comp[j]->getA())[1];
        }
        for (int i = 1; i <= k; i++) {
          if (comp[j]->get_type() == MRG && i <= comp[j]->getK()) A[i] += a[i-1]*n*modulo/mod_gen;
          if (comp[j]->get_type() == MWC && i <= 1) A[i] += a[i-1]*n*modulo/mod_gen;
        }
      }
      // Modulo only once
      for (int i = 1; i <= k; i++) {
        A[i] %= modulo;
      }

      return new MRGLattice<Int, Dbl>(modulo, A, maxDim, k, FULL);
    }

  extern template MRGLattice<std::int64_t, double>* getLatCombo(std::vector<MRGComponent<std::int64_t>*>& comp, int maxDim);
  extern template MRGLattice<NTL::ZZ, double>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>*>& comp, int maxDim);
  extern template MRGLattice<NTL::ZZ, NTL::RR>* getLatCombo(std::vector<MRGComponent<NTL::ZZ>*>& comp, int maxDim);

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
         * Ideally, `lat` has been initialized with `getLatCombo`.
         * */
        ComboLattice(std::vector<MRGComponent<Int>*>& comp,
            MRGLattice<Int, Dbl>& lat);

        /**
         * Copy constructor
         * */
        ComboLattice(const ComboLattice<Int, Dbl>& lat);

        /**
         * Destructor.
         * Deallocate memory to `m_comp`.
         * */
        ~ComboLattice();

        /**
         * Prints a string describing all the componenets and the equivalent MRG.
         * */
        std::string toString() const override;

      private:
        std::vector<MRGComponent<Int>*> m_comp;

        /**
         * The number of components stored in this object.
         * */
        int m_number;
    }; // end class ComboLattice

  //============================================================================

  template<typename Int, typename Dbl>
    ComboLattice<Int, Dbl>::ComboLattice(std::vector<MRGComponent<Int>*>& comp,
        MRGLattice<Int, Dbl>& lat) : MRGLattice<Int, Dbl>(lat){
      for (unsigned int i = 0; i < comp.size(); i++) {
        m_comp.push_back(new MRGComponent<Int>(*comp[i]));
      }
      m_number = comp.size();
    }

  //============================================================================

  template<typename Int, typename Dbl>
    ComboLattice<Int, Dbl>::ComboLattice(const ComboLattice<Int, Dbl>& lat) :
      MRGLattice<Int, Dbl>(lat){
      for (unsigned int i = 0; i < lat.m_comp.size(); i++) {
        m_comp.push_back(new MRGComponent<Int>(*lat.m_comp[i]));
      }
      m_number = lat.m_number;
    }

  //============================================================================

  template<typename Int, typename Dbl>
    ComboLattice<Int, Dbl>::~ComboLattice() {
      for (unsigned int i = 0; i < m_comp.size(); i++) {
        delete m_comp[i];
      }
    }

  //============================================================================

  template<typename Int, typename Dbl>
    std::string ComboLattice<Int, Dbl>::toString() const {
      std::ostringstream out;
      for (int i = 0; i < m_number; i++) {
        Int m;
        if (m_comp[i]->get_type() == MRG) m = m_comp[i]->getM();
        if (m_comp[i]->get_type() == MWC) m = m_comp[i]->m_MWCb;
        out << "Component " << i+1 << "\nm = " << m << "\nk = "
          << m_comp[i]->getK() << "\na = " << m_comp[i]->getA() << "\n\n";
      }
      out << "Equivalent to\n";
      out << "m = " << this->m_modulo << "\nk = " << this->m_order << "\n";
      for (int i = 0; i < this->m_order; i++) {
        out << "a_" << i+1 << " = " << this->m_aCoef[i];
        out << "\n";
      }
      return out.str();
    }

  //============================================================================

  extern template class ComboLattice<std::int64_t, double>;
  extern template class ComboLattice<NTL::ZZ, double>;
  extern template class ComboLattice<NTL::ZZ, NTL::RR>;

} // end namespace LatMRG
#endif

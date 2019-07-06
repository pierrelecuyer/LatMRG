#ifndef LATMRG_LCGCARRYLATTICE_H
#define LATMRG_LCGCARRYLATTICE_H

#include "latmrg/MRGLattice.h"
#include "latmrg/IntFactorization.h"

namespace LatMRG {

  /**
   * This class represents the lattice a of LCG with a carry defined by the
   * recurrence
   * \f[
   *   x_n = a x_{n-1} +c \pmod m.
   * \f]
   * The building
   * strategy for the lattice is technically the same, but MRGLattice does not
   * feature the possibility to add one becausae it is seldom used for a MRG.
   * The class overrides the basis construction methods and implements a method
   * to check if the parameter choice \f$(a,c,m)\f$ defines a full period LCG
   * with carry.
   */
  template<typename Int, typename Dbl>
    class LCGCarryLattice: public MRGLattice<Int, Dbl> {
      public:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;

        /**
         * Initializes a LCGCarry lattice for the recurrence
         * \f[
         *   x_n = a x_{n-1} +c \pmod m.
         * \f]
         * */
        LCGCarryLattice (const Int & a, const Int & c, const Int & m);

        /**
         * Copy constructor.
         * */
        LCGCarryLattice (const LCGCarryLattice<Int, Dbl> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        LCGCarryLattice<Int, Dbl> & operator= (const LCGCarryLattice<Int, Dbl> & Lat);
                
        /**
         * Gets the parameters of the lattice `(a, c, m)` in a formated string.
         */
        std::string toString() const override;

        /**
         * Builds the lattice basis for dimension `d`.
         * */
        void buildBasis(int d) override;

        /**
         * Increases the dimension of the lattice by one.
         * */
        void incDimBasis() override;

        /**
         *
         * */

        /**
         * This is a basic method to check if the MWC generator described by
         * `b` and `e` has full period.
         *
         * \todo finish this
         * */
        static bool fullPeriod(const Int& a, const Int& c, const Int& m) {
          if (NTL::GCD(c, m) > 1) return false;
          IntFactorization<Int> m_factor(m);
          auto fact_list = m_factor.getFactorList();
          for (auto iter = fact_list.begin(); iter != fact_list.end(); iter++) {
            if (((a - Int(1)) % (*iter).getFactor()) != Int(0)) return false;
          }
          if ((m % Int(4)) != Int(0)) return (((a-Int(1))%Int(4)) == 0);
          return true;
        }

      private:

        /**
         * Storing the carry of the recurrence.
         * */
        Int m_carry;
    }; // End class declaration

  //===========================================================================

  template<typename Int, typename Dbl>
    LCGCarryLattice<Int, Dbl>::LCGCarryLattice (const Int & a, const Int & c, const Int & m, int maxDim):
      MRGLattice(m, a, maxDim, FULL) {
        m_carry = c;
      }

  //===========================================================================

  template<typename Int, typename Dbl>
    LCGCarryLattice<Int, Dbl>::LCGCarryLattice (const LCGCarryLattice<Int, Dbl> & Lat):
      MRGLattice(Lat) {
        m_carry = Lat.m_carry;
      }

  //===========================================================================

  template<typename Int, typename Dbl>
    LCGCarryLattice<Int, Dbl> & LCGCarryLattice<Int, Dbl>::operator= (
        const LCGCarryLattice<Int, Dbl> & Lat) {
      MRGLattice<Int, Dbl>::operator=(Lat);
      m_carry = Lat.m_carry;
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    std::string LCGCarryLattice<Int, Dbl>::toString() const {
      std::ostringstream out;
      out << "a = " << m_aCoef[0] << "\n";
      out << "c = " << m_carry << "\n";
      return out.str ();
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    void LCGCarryLattice<Int, Dbl>::buildBasis(int d) {
      this->m_basis.SetDims(1,1);
      this->m_dualbasis.SetDims(1,1);
      this->setDim(1);

      this->m_basis[0][0] = Int(1);
      this->m_dualbasis[0][0] = this->m_modulo;

      for (int i = 0; i<d-1; i++)
        incDimBasis();

  //===========================================================================

  template<typename Int, typename Dbl>
    void LCGCarryLattice<Int, Dbl>::incDimBasis() {

      LatticeTester::IntLattice<Int, Int, Dbl, Dbl>::incDim();
      const int dim = this->getDim();

      IntMat temp1(this->m_basis);
      IntMat temp2(this->m_dualbasis);

      for (int i = 0; i < dim; i++) {
        this->m_basis[i][dim-1] = 0;
        this->m_basis[dim-1][i] = 0;
        this->m_dualbasis[i][dim-1] = 0;
        this->m_dualbasis[dim-1][i] = 0;
      }
      this->m_basis[0][dim-1] = (m_aCoef[0]*this->m_basis[0][dim-2] + m_carry) % this->m_modulo;
      this->m_dualbasis[dim-1][0] = -this->m_basis[0][dim-1];
      this->m_basis[dim-1][dim-1] = this->m_modulo;
      this->m_dualbasis[dim-1][dim-1] = 1;

      this->setNegativeNorm();
      this->setDualNegativeNorm();
    }

  //============================================================================
  //The types combinations supported in the library

  extern template class LCGCarryLattice<std::int64_t, double>;
  extern template class LCGCarryLattice<NTL::ZZ, double>;
  extern template class LCGCarryLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif


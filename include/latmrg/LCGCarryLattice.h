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
         * Destructor.
         */
        ~LCGCarryLattice();

        /**
         * Cleans and releases memory used by this object.
         */
        void kill();

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
        void incDim() override;

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

  //============================================================================
  //The types combinations supported in the library

  extern template class LCGCarryLattice<std::int64_t, double>;
  extern template class LCGCarryLattice<NTL::ZZ, double>;
  extern template class LCGCarryLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif


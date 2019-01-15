#ifndef MWCLATTICE_H
#define MWCLATTICE_H

#include "latticetester/Const.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLattice.h"

#include "latmrg/Const.h"
#include "latmrg/MRGComponent.h"
#include "latmrg/MRGLattice.h"

#include <string>

namespace LatMRG {

  /**
   * This class represents the lattice associated to a Multiply-with-carry (MWC)
   * random number generator. A MWC generator is defined by a recurrence of the
   * form
   * \f{align}{
   *    x_n & = (e_1 x_{n-1} + \cdots + e_k x_{n-k} + c_{n-1})d\ \mathrm{mod} \ b \\
   *    c_n & = \lfloor (e_0 x_n + e_1 x_{n-1} + \cdots + e_k x_{n-k} + c_{n-1} )/b \rfloor \\
   *    u_n & = \sum_{i=1}^\infty x_{n+i-1} b^{-i}
   * \f}
   * This generator can then be reprensented as an LCG and this class simply
   * builds the correct MRGLattice object corresponding to a specific MWC
   * generator. All the functions and attributes inherited by this class
   * work as they would be expected to work on the lattice corresponding to that
   * MRGLattice. For example, this means that the vector m_aCoef stores the
   * coefficient of the LCG instead of the coefficients of the MWC generator.
   *
   * This class simply implements a constructor and and the functions to compute
   * the LCG equivalent to the MRG.
   */
  template<typename Int, typename Dbl>
    class MWCLattice: public MRGLattice<Int, Dbl> {
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;
      public:

        /**
         * b is the modulo.
         * e is the vector of coefficients.
         * k is the order of the recurrence.
         *
         * It is recommended to verify that the parameters passed to this
         * constructor work before using it, because the program will crash of
         * `validate(b, e) != 0`.
         */
        MWCLattice (const Int & b, const IntVec & e, int k);

        /**
         * Another constructor with just `b` and `m`. Of course, `b` and `m`
         * have to be valid.
         * */
        MWCLattice (const Int & b, const Int & m);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        MWCLattice (const MWCLattice<Int, Dbl> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        MWCLattice<Int, Dbl> & operator= (const MWCLattice<Int, Dbl> & Lat);

        /**
         * Destructor.
         */
        ~MWCLattice();

        /**
         * Cleans and releases memory used by this object.
         */
        void kill();

        /**
         * Reimplement this.
         */
        std::string toStringCoef() const;

        /**
         * The modulo of the MWC generator this object represents.
         * */
        const Int& getMWCmod() const {return this->m_MWCmod;}

        /**
         * The order of the MWC generator. An order of `k` means that the
         * generator has `k+1` coefficients.
         * */
        int getMWCorder() const {return this->m_MWCorder;}

        /**
         * The coefficients of the MWC generator. We have that
         * `m_eCoef[i] = e_i`. The MWC equations mean that if the order is `k`,
         * the generator has `k+1` coefficients.
         * */
        const IntVec& geteCoef() const {return this->m_eCoef;}

        /**
         * This checks if `b` and `e` are suitable parameters for a MWC
         * generator.
         * It is recommended to call this function before building a MWCLattice
         * object because the program will exit if the condition verified here
         * is not met.
         *
         * This basically just checks that `gcd(b, e[0]) = 1` and returns `1` if
         * it is `false` and `0` if it checks out.
         * */
        static int validate(const Int& b, const IntVec& e) {
          if(NTL::GCD(b, e[0]) != 1) return 1;
          return 0;
        }

        /**
         * This is a basic method to check if the MWC generator described by
         * `b` and `e` has full period.
         * */
        static int fullPeriod(const Int& b, const IntVec& e) {
          if(MWCLattice<Int, Dbl>::validate(b, e)) return 1;
          return 0;
        }

      private:

        /**
         * Returns the coefficient for an MWC with coefficients in `e` and
         * modulo `b`.
         * */
        IntVec LCGCoeff(const Int& b, const IntVec& e);

        /**
         * Returns the modulo for an MWC with coefficients in `e` and
         * modulo `b`.
         * */
        Int LCGMod(const Int& b, const IntVec& e);

        /**
         * The modulo of the MWC generator this object represents.
         * */
        Int m_MWCmod;

        /**
         * The order of the MWC generator. An order of `k` means that the
         * generator has `k+1` coefficients.
         * */
        int m_MWCorder;

        /**
         * The coefficients of the MWC generator. We have that
         * `m_eCoef[i] = e_i`. The MWC equations mean that if the order is `k`,
         * the generator has `k+1` coefficients.
         * */
        IntVec m_eCoef;

    }; // End class declaration

  //===========================================================================

  template<typename Int, typename Dbl>
    MWCLattice<Int, Dbl>::MWCLattice(const Int & b, const IntVec & e, int k):
      MRGLattice<Int, Dbl>(this->LCGMod(b, e), this->LCGCoeff(b,e), 1, 1, FULL)
  {
    m_MWCmod = b;
    m_MWCorder = k;
    m_eCoef.SetLength(k+1);
    for (int i = 0; i < k+1; i++)
      m_eCoef[i] = e[i];
  }

  //===========================================================================

  template<typename Int, typename Dbl>
    MWCLattice<Int, Dbl>::MWCLattice(const Int & b, const Int & m):
      MRGLattice<Int, Dbl>(m, NTL::InvMod(b, m), 1, FULL)
  {
    m_MWCmod = b;
    m_MWCorder = 0;
    // Even though this is not needed, this constructor should set the
    // coefficients of the MWC generator to the representation of m in base b.
    m_eCoef.SetLength(1);
  }

  //===========================================================================

  template<typename Int, typename Dbl>
    MWCLattice<Int, Dbl>::MWCLattice(const MWCLattice<Int, Dbl> &lat):
      MRGLattice<Int, Dbl> (lat)
  {
    m_MWCmod = lat.getMWCmod();
    m_MWCorder = lat.getMWCorder();
    m_eCoef.SetLength(m_MWCorder+1);
    for (int i = 0; i < m_MWCorder+1; i++)
      m_eCoef[i] = lat.geteCoef()[i];
  }


  //===========================================================================

  template<typename Int, typename Dbl>
    MWCLattice<Int, Dbl> & MWCLattice<Int, Dbl>::operator= (const MWCLattice<Int, Dbl> & lat)
    {
      if (this == &lat)
        return *this;
      (MRGLattice<Int, Dbl>) *this = (MRGLattice<Int, Dbl>)lat;
      m_MWCmod = lat.getMWCmod();
      m_MWCorder = lat.getMWCorder();
      m_eCoef.SetLength(m_MWCorder+1);
      for (int i = 0; i < m_MWCorder+1; i++)
        m_eCoef[i] = lat.geteCoef()[i];
      return *this;
    }

  //============================================================================

  template<typename Int, typename Dbl>
    MWCLattice<Int, Dbl>::~MWCLattice ()
    {
      kill();
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void MWCLattice<Int, Dbl>::kill()
    {
      MRGLattice<Int, Dbl>::kill();
      m_eCoef.kill();
    }

  //============================================================================

  template<typename Int, typename Dbl>
    std::string MWCLattice<Int, Dbl>::toStringCoef () const
    {
      std::ostringstream out;
      out << "LCG coefficient: ";
      //out << "[ ";
      for (int i = 0; i < this->m_order; i++)
        out << this->m_aCoef[i] << "  ";
      //out << "]";
      out << "MWC coefficients: ";
      for (int i = 0; i <= this->m_MWCorder; i++)
        out << m_eCoef[i] << "  ";
      return out.str ();
    }

  //============================================================================

  template<typename Int, typename Dbl>
    typename MWCLattice<Int,Dbl>::IntVec MWCLattice<Int, Dbl>::LCGCoeff(const Int& b, const IntVec& e){
      Int mult = this->LCGMod(b,e);
      Int a = NTL::InvMod(b, mult);
      IntVec coeff;
      coeff.SetLength(2);
      coeff[1] = a;
      return coeff;
    }

  //============================================================================

  template<typename Int, typename Dbl>
    Int MWCLattice<Int, Dbl>::LCGMod(const Int& b, const IntVec& e){
      Int m(0);
      for(int i = 0; i <= e.length(); i++) {
        m += e[i] * NTL::power(b, i);
      }
      return m;
    }

  //============================================================================
  //The types combinations supported in the library

  extern template class MWCLattice<std::int64_t, double>;
  extern template class MWCLattice<NTL::ZZ, double>;
  extern template class MWCLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif


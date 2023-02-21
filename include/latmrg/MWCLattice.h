#ifndef LATMRG_MWCLATTICE_H
#define LATMRG_MWCLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"

#include "latmrg/EnumTypes.h"
#include "latmrg/MRGPeriod.h"
#include "latmrg/MRGLattice.h"

#include <string>

/**
 * Small functions to give the modulo and the coefficient of the LCG generator
 * equivalent to a MWC generator with modulo b and coefficients e.
 * */
namespace MWCEquiv {
  /**
   * Returns the modulo for an MWC with coefficients in `e` and
   * modulo `b`.
   * */
  template<typename Int>
    Int LCGMod(const Int& b, const NTL::vector<Int>& e){
      Int m(0);
      for(int i = 0; i < e.length(); i++) {
        m += e[i] * NTL::power(b, i);
      }
      return m;
    }

  //============================================================================

  /**
   * Returns the coefficient for an MWC with coefficients in `e` and
   * modulo `b`.
   * */
  template<typename Int>
    NTL::vector<Int> LCGCoeff(const Int& b, const NTL::vector<Int>& e){
      Int mult = LCGMod(b,e);
      Int a = NTL::InvMod(b, mult);
      NTL::vector<Int> coeff;
      coeff.SetLength(2);
      coeff[1] = a;
      return coeff;
    }
}

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
  template<typename Int, typename Real>
    class MWCLattice: public MRGLattice<Int, Real> {
      public:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;

        /**
         * b is the modulo.
         * e is the vector of coefficients.
         * k is the order of the recurrence.
         *
         * It is recommended to verify that the parameters passed to this
         * constructor work before using it, because the program will crash of
         * `validate(b, e) != 0`.
         */
        MWCLattice (const Int & b, const IntVec & e, int k, int maxDim = 1);

        /**
         * Another constructor with just `b` and `m`. Of course, `b` and `m`
         * have to be valid.
         * */
        MWCLattice (const Int & b, const Int & m, int maxDim = 1);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        MWCLattice (const MWCLattice<Int, Real> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        MWCLattice<Int, Real> & operator= (const MWCLattice<Int, Real> & Lat);

        /**
         * Destructor.
         */
        ~MWCLattice();

        /**
         * Cleans and releases memory used by this object.
         */
        void kill();

        /**
         * Gets the coefficients of the MRG that spawns the lattice in a string.
         */
        std::string toStringCoef() const;

        /**
         * Gets all the information on the lattice in a string. This contains
         * the type of generator and the coefficients.
         */
        std::string toString() const override;

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
         *
         * \todo finish this
         * */
        static int fullPeriod(const Int& b, const IntVec& e) {
          if(MWCLattice<Int, Real>::validate(b, e)) return 1;
          return 0;
        }

      private:

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

  template<typename Int, typename Real>
    MWCLattice<Int, Real>::MWCLattice(const Int & b, const IntVec & e, int k, int maxDim):
      MRGLattice<Int, Real>(MWCEquiv::LCGMod(b, e), MWCEquiv::LCGCoeff(b,e), maxDim, 1, FULL)
  {
    m_MWCmod = b;
    m_MWCorder = k;
    m_eCoef.SetLength(k+1);
    for (int i = 0; i < k+1; i++)
      m_eCoef[i] = e[i];
  }

  //===========================================================================

  template<typename Int, typename Real>
    MWCLattice<Int, Real>::MWCLattice(const Int & b, const Int & m, int maxDim):
      MRGLattice<Int, Real>(m, NTL::InvMod(b, m), maxDim, FULL)
  {
    m_MWCmod = Int(b);
    // This is not needed in reality, but this tries to compute the coefficients
    // of the MWC generator. This is imperfect because we do not know if the
    // coefficients are negative.
    Int modulo(m);
    int k = -1;
    m_eCoef.SetLength(0);
    while (modulo != 0) {
      k++;
      Int rest = modulo % b;
      if (rest > b/2) rest = rest-b;
      m_eCoef.append(rest);
      modulo -= rest;
      modulo /= b;
    }
    m_MWCorder = k;
  }

  //===========================================================================

  template<typename Int, typename Real>
    MWCLattice<Int, Real>::MWCLattice(const MWCLattice<Int, Real> &lat):
      MRGLattice<Int, Real> (lat)
  {
    m_MWCmod = Int(lat.getMWCmod());
    m_MWCorder = lat.getMWCorder();
    m_eCoef.SetLength(m_MWCorder+1);
    for (int i = 0; i < m_MWCorder+1; i++)
      m_eCoef[i] = Int(lat.geteCoef()[i]);
  }


  //===========================================================================

  template<typename Int, typename Real>
    MWCLattice<Int, Real> & MWCLattice<Int, Real>::operator= (const MWCLattice<Int, Real> & lat)
    {
      if (this == &lat)
        return *this;
      (MRGLattice<Int, Real>) *this = (MRGLattice<Int, Real>)lat;
      m_MWCmod = lat.getMWCmod();
      m_MWCorder = lat.getMWCorder();
      m_eCoef.SetLength(m_MWCorder+1);
      for (int i = 0; i < m_MWCorder+1; i++)
        m_eCoef[i] = lat.geteCoef()[i];
      return *this;
    }

  //============================================================================

  template<typename Int, typename Real>
    MWCLattice<Int, Real>::~MWCLattice ()
    {
      kill();
    }


  //===========================================================================

  template<typename Int, typename Real>
    void MWCLattice<Int, Real>::kill()
    {
      MRGLattice<Int, Real>::kill();
      m_eCoef.kill();
    }

  //============================================================================

  template<typename Int, typename Real>
    std::string MWCLattice<Int, Real>::toStringCoef () const
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

  template<typename Int, typename Real>
    std::string MWCLattice<Int, Real>::toString () const
    {
      std::ostringstream out;
      for (int i = 0; i <= this->m_MWCorder; i++) {
        out << "a_" << i << " = " << m_eCoef[i];
        out << "\n";
      }

      out << "\nLCG equivalent:\n";
      out << "m = " << this->m_modulo << "\n";
      out << "a = " << this->m_aCoef[0] << "\n";
      return out.str ();
    }

  //============================================================================
  //The types combinations supported in the library

  extern template class MWCLattice<std::int64_t, double>;
  extern template class MWCLattice<NTL::ZZ, double>;
  extern template class MWCLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif


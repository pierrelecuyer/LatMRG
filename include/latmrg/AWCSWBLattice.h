#ifndef LATMRG_AWCSWBLATTICE_H
#define LATMRG_AWCSWBLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"

#include "latmrg/EnumTypes.h"
#include "latmrg/MRGPeriod.h"
#include "latmrg/MRGLattice.h"

#include <string>

namespace {
  /**
   * Returns the modulo for the LCG associated with the AWC/SWB describeb by the
   * parameters.
   * */
  template<typename Int>
    Int LCGMod(const Int& b, int r, int s, int mode){
      Int m(0);
      if (mode == 0) {
        m = NTL::power(b, r) + NTL::power(b, s) - 1;
      } else if (mode == 1) {
        m = NTL::power(b, r) + NTL::power(b, s) + 1;
      } else if (mode == 2) {
        m = NTL::power(b, r) - NTL::power(b, s) + 1;
      } else if (mode == 3) {
        m = NTL::power(b, r) - NTL::power(b, s) - 1;
      }
      return m;
    }

  //============================================================================

  /**
   * Returns the coefficient for the LCG associated with the AWC/SWB describeb
   * by the parameters.
   * */
  template<typename Int>
    NTL::vector<Int> LCGCoeff(const Int& b, int r, int s, int mode){
      Int mult = LCGMod(b, r, s, mode);
      Int a = NTL::InvMod(b, mult);
      NTL::vector<Int> coeff;
      coeff.SetLength(2);
      coeff[1] = a;
      return coeff;
    }
}

namespace LatMRG {

  /**
   * This class represents the lattice associated with either an AWC or a SWB
   * random number generator. This class is based on the article \cite rTEZ93a
   * and uses the exact same representation, except that each case is
   * represented here by an integer ranging from 0 to 3. It is important to note
   * that when interracting with other classes of LatMRG, this will most likely
   * act as a simple MRGLattice and only the LCG equivalent to the AWC/SWB will
   * be considered. When using this class, `mode = 0` means that the generator
   * is of the form AWC, `mode = 1` is equivalent to AWC-c, `mode = 2` is SWB-I
   * and `mode = 3` is SWB-II.
   */
  template<typename Int, typename Real>
    class AWCSWBLattice: public MRGLattice<Int, Real> {
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;
      public:

        /**
         * b is the modulo.
         * We expect that `r > s` if this is not the case, then `r` and `s` are
         * reversed.
         */
        AWCSWBLattice (const Int & b, int r, int s, int mode);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        AWCSWBLattice (const AWCSWBLattice<Int, Real> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        AWCSWBLattice<Int, Real> & operator= (const AWCSWBLattice<Int, Real> & Lat);

        /**
         * Destructor.
         */
        ~AWCSWBLattice();

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
        const Int& getAWCSWBmod() const {return this->m_bModulo;}

        /**
         * The order of the MWC generator. An order of `k` means that the
         * generator has `k+1` coefficients.
         * */
        int getAWCSWBorder() const {return this->m_subOrder;}

      private:

        /**
         * The modulo of the AWC/SWB generator this object represents.
         * */
        Int m_bModulo;

        /**
         * The `r` parameter to this type of generator.
         * */
        int m_subOrder;

        /**
         * Parameter to this AWC/SWB generator.
         * */
        int m_s;

        /**
         * The type of generator this is.
         * */
        int m_mode;

    }; // End class declaration

  //===========================================================================

  template<typename Int, typename Real>
    AWCSWBLattice<Int, Real>::AWCSWBLattice(const Int & b, int r, int s, int mode):
      MRGLattice<Int, Real>(LCGMod(b,r,s,mode), LCGCoeff(b,r,s,mode), 1, 1, FULL)
  {
    if (r > s) {
      m_subOrder = r;
      m_s = s;
    } else {
      m_subOrder = s;
      m_s = r;
    }
    m_bModulo = b;
    m_mode = mode;
  }

  //===========================================================================

  template<typename Int, typename Real>
    AWCSWBLattice<Int, Real>::AWCSWBLattice(const AWCSWBLattice<Int, Real> &lat):
      MRGLattice<Int, Real> (lat)
  {
    m_bModulo = Int(lat.getAWCSWBmod());
    m_subOrder = lat.getAWCSWBorder();
    m_s = lat.m_s;
    m_mode = lat.m_mode;
  }


  //===========================================================================

  template<typename Int, typename Real>
    AWCSWBLattice<Int, Real> & AWCSWBLattice<Int, Real>::operator= (const AWCSWBLattice<Int, Real> & lat)
    {
      if (this == &lat)
        return *this;
      (MRGLattice<Int, Real>) *this = (MRGLattice<Int, Real>)lat;
      m_bModulo = lat.getAWCSWBmod();
      m_subOrder = lat.getAWCSWBorder();
    m_s = lat.m_s;
    m_mode = lat.m_mode;
      return *this;
    }

  //============================================================================

  template<typename Int, typename Real>
    AWCSWBLattice<Int, Real>::~AWCSWBLattice ()
    {
      kill();
    }


  //===========================================================================

  template<typename Int, typename Real>
    void AWCSWBLattice<Int, Real>::kill()
    {
      MRGLattice<Int, Real>::kill();
    }

  //============================================================================

  template<typename Int, typename Real>
    std::string AWCSWBLattice<Int, Real>::toStringCoef () const
    {
      std::ostringstream out;
      out << "LCG coefficient: ";
      //out << "[ ";
      for (int i = 0; i < this->m_order; i++)
        out << this->m_aCoef[i] << "  ";
      //out << "]";
      out << "AWC/SWB info:\n";
      out << "r: " << this->m_subOrder << "\ns: " << this->m_s << "\nmode: " <<
        this->m_mode;
      return out.str ();
    }

  //============================================================================


  //============================================================================
  //The types combinations supported in the library

  extern template class AWCSWBLattice<std::int64_t, double>;
  extern template class AWCSWBLattice<NTL::ZZ, double>;
  extern template class AWCSWBLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif


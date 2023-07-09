#ifndef LATMRG_LCGLATTICE_H
#define LATMRG_LCGLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"


namespace LatMRG {

  /**
   * This class implements lattice bases built from a Korobov lattice rule. For
   * a given \f$a\f$, a Korobov lattice basis is formed as follows:
   * \f[
   * \mathbf{b_1} = (1, a, a^2, …, a^{d-1}),\quad
   * \mathbf{b_2} = (0, n, 0, …, 0),\quad…,\quad
   * \mathbf{b_d} = (0, …, 0, n).
   * \f]
   *
   * \remark **Pierre:** Reprogrammer \c incDim de façon efficace comme dans
   * \c MRGLattice
   *
   */
template<typename Int, typename Real>
    class LCGLattice: public LatticeTester::IntLatticeExt<Int, Real> {
      private:
        typedef NTL::matrix<Int> IntMat;
      public:

        /**
         * Constructs a Korobov lattice with \f$n\f$ points, maximal dimension
         * `maxDim` using the norm `norm`.
         */
        LCGLattice (const Int & n, const Int & a, int maxDim,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Constructor. Same as above, except the lattice is formed as follow:
         * \f[
         * \mathbf{b_1} = (a^t, a^{t+1}, a^{t+2}, …, a^{t+d-1}),\qquad
         * \mathbf{b_2} = (0, n, 0, …, 0),\qquad…,\qquad
         * \mathbf{b_d} = (0, …, 0, n).
         * \f]
         */
        LCGLattice (const Int & n, const Int & a, int dim, int t,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Copy constructor.
         */
        LCGLattice (const LCGLattice & Lat);

        /**
         * Assigns `Lat` to this object.
         */
        LCGLattice & operator= (const LCGLattice & Lat);

        /**
         * Destructor.
         */
        ~LCGLattice();

        /**
         * Returns the multiplier \f$a\f$ as a string.
         */
        std::string toStringCoef() const;

        /**
         * Builds the basis in dimension \f$d\f$.
         */
        void buildBasis (int d);

        /**
         * Increments the dimension of the basis by 1.
         */
        void incDim();

        /**
         * Increments the dimension of the basis by 1 by rebuilding the basis
         * from scratch. This is very slow. It can be used for verification of
         * the fast `incDim` method above.
         */
        void incDimSlow();

      protected:

        /**
         * The multiplier of the Korobov lattice rule.
         */
        Int m_a;

        /**
         * The shift applied to the lattice rule.
         */
        int m_shift;

        /**
         * Initialization.
         */
        void init();
    }; // End class declaration

  template<typename Int, typename Real>
    LCGLattice<Int, Real>::LCGLattice (const Int & n, const Int & a,
        int maxDim, LatticeTester::NormType norm) :
      LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(n, 0, maxDim, norm)
  {
    m_a = a;
    m_shift = 0;
    init();
  }


  //===========================================================================

  template<typename Int, typename Real>
    LCGLattice<Int, Real>::LCGLattice (const Int & n, const Int & a,
        int maxDim, int t, LatticeTester::NormType norm) :
      LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(n, 0, maxDim, norm)
  {
    m_a = a;
    m_shift = t;
    init();
  }


  //===========================================================================

  template<typename Int, typename Real>
    LCGLattice<Int, Real>::~LCGLattice()
    {
    }


  //=========================================================================

  template<typename Int, typename Real>
    LCGLattice<Int, Real>::LCGLattice (const LCGLattice & lat):
      LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(
    		  lat.m_modulo, 0, lat.getDim (), lat.getNorm ())
  {
    m_a = lat.m_a;
    m_shift = lat.m_shift;
  }


  //=========================================================================

  template<typename Int, typename Real>
    LCGLattice<Int, Real> & LCGLattice<Int, Real>::operator=
    (const LCGLattice & lat)
    {
      if (this == &lat)
        return *this;
      this->m_dim = lat.m_dim;
      this->copyBasis(lat);
      this->m_order = lat.m_order;
      m_a = lat.m_a;
      m_shift = lat.m_shift;
      return *this;
    }


  //===========================================================================

  template<typename Int, typename Real>
    void LCGLattice<Int, Real>::init()
    {
      //Erwan   IntLattice::init();
      //   double temp;
      //   conv (temp, m_m);
      //   m_lgVolDual2[1] = 2.0 * Lg(temp);

      //Erwan   for (int r = m_order + 1; r <= getMaxDim(); r++)
      //Erwan      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
    }


  //===========================================================================

  template<typename Int, typename Real>
    std::string LCGLattice<Int, Real>::toStringCoef () const
    {
      std::ostringstream out;
      out << m_a;
      return out.str ();
    }


  //===========================================================================

  template<typename Int, typename Real>
    void LCGLattice<Int, Real>::buildBasis(int d)
    {
      //assert(d <= getMaxDim());
      this->setDim(d);
      Int tmp;
      // the .0 double-literal is needed on some systems to prevent compilation error
      NTL::conv(tmp, 1.0);

      for (int i = 0; i < m_shift; i++) {
        tmp *= m_a;
        tmp %= this->m_modulo;
      }

      for (int j = 1; j <= d; j++) {
        // V[1][j] = tmp % m;
        this->m_basis (0, j) = tmp;
        tmp *= m_a;
        tmp %= this->m_modulo;
      }

      for (int i = 2; i <= d; i++) {
        for (int j = 1; j <= d; j++) {
          if (i == j)
            this->m_basis (i, j) = this->m_modulo;
          else
            this->m_basis (i, j) = 0;
        }
      }

      // Build dual basis
      LatticeTester::CalcDual<IntMat>(this->m_basis, this->m_dualbasis, d, this->m_modulo);
    }


  //===========================================================================

  template<typename Int, typename Real>
    void LCGLattice<Int, Real>::incDimSlow()
    {
      // Temporaire: très lent. Reprogrammer.
      int d = this->getDim();
      buildBasis(d + 1);
      this->setNegativeNorm();
      this->setDualNegativeNorm();
    }


  //===========================================================================

  template<typename Int, typename Real>
    void LCGLattice<Int, Real>::incDim()
    {
      Int tmp1, tmp2, tmp3;
      NTL::vector<Int> vectmp1;// working variables
      LatticeTester::IntLatticeExt<Int, Real>::incDim(); // Increments the dimenson by 1
      const int dim = this->getDim();  // New dimension

      vectmp1.resize(dim);
      for (int i = 1; i < dim-1; i++) {
        NTL::conv (tmp2, this->m_basis (i, dim - 2));
        tmp1 = tmp2 * m_a;
        LatticeTester::Modulo (tmp1, this->m_modulo, tmp1);
        this->m_basis (i, dim) = vectmp1(i) =  tmp1;    // Erwan m_vSI (0, i) = tmp1;
      }

      NTL::matrix_row<IntMat> row1(this->m_basis, dim - 2);
      LatticeTester::SetZero (row1, dim - 2);

      for (int i = 0; i < dim-1; i++)
        this->m_basis (dim-1, i) = 0;
      this->m_basis (dim-1, dim-1) = this->m_modulo;

      for (int i = 0; i < dim-1; i++)
        this->m_dualbasis (i, dim-1) = 0;
      this->m_dualbasis (dim-1, dim-1) = 1;

      for (int j = 1; j < dim; j++) {
        NTL::clear (tmp1);
        for (int i = 1; i < dim; i++) {
          tmp2 = this->m_dualbasis (i, j);
          tmp2 *= vectmp1 (i);
          tmp1 -= tmp2;
        }
        LatticeTester::Quotient (tmp1, this->m_modulo, tmp3);
        this->m_dualbasis (dim, j) = tmp3;
      }

      this->setNegativeNorm();
      this->setDualNegativeNorm();
    }

  extern template class LCGLattice<std::int64_t, double>;
  extern template class LCGLattice<NTL::ZZ, double>;
  extern template class LCGLattice<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif

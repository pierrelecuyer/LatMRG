#ifndef MRGLATTICE_H
#define MRGLATTICE_H
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLattice.h"

#include "latmrg/Const.h"
#include "latmrg/MRGComponent.h"

#include <string>


namespace LatMRG {

  /**
   * This class implements lattice basis built from multiple recursive linear
   * congruential generators (MRGs). One must first call the constructor with a
   * given congruence modulus \f$m\f$, a given order \f$k\f$ for the
   * recurrence, and a maximal dimension for the basis. One must then build the
   * lattice basis associated to a vector of multipliers for a given dimension.
   * Each MRG is defined by a vector of multipliers \f$A\f$, where \f$A[i]\f$
   * represents \f$a_i\f$. This MRG satisfies the recurrence
   * \f[
   *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mod m.
   * \f]
   */
  template<typename Int>
    class MRGLattice: public LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal> {
      public:

        /**
         * Constructor with modulus of congruence \f$m\f$, order of the 
         * recurrence \f$k\f$, multipliers \f$a\f$, maximal dimension `MaxDim`, 
         * and lattice type `Latt`. Vectors and (square) matrices of the basis 
         * have maximal dimension `maxDim`, and the indices of vectors and 
         * matrices vary from dimension 1 to `maxDim`. The norm to be used for 
         * the basis vectors is `norm`.
         */
        MRGLattice (const Int & m, const MVect & a, int maxDim, int k, 
            LatticeType latt, 
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * As in the constructor above but the basis is built for the lacunary
         * indices `lac`.
         */
        MRGLattice (const Int & m, const MVect & a, int maxDim, int k, 
            BVect & lac, LatticeType latt, 
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        MRGLattice (const MRGLattice<Int> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        MRGLattice<Int> & operator= (const MRGLattice<Int> & Lat);

        /**
         * Destructor.
         */
        ~MRGLattice();

        /**
         * Cleans and releases memory used by this object.
         */
        void kill();

        /**
         * Builds the basis in dimension \f$d\f$.
         */
        virtual void buildBasis (int d);

        /**
         * Increments the dimension of the basis by 1 by calling either
         * `incDimBasis` or `incDimLaBasis`.
         */
        virtual void incDim();

        /**
         * Returns `true` for the case of lacunary indices, returns `false` for
         * non-lacunary indices.
         */
        bool isLacunary() const { return m_lacunaryFlag; }

        /**
         * Returns the \f$j\f$-th lacunary index.
         */
        BScal & getLac (int j);

        /**
         * Sets the lacunary indices for this lattice to `lat`.
         */
        virtual void setLac (const LatticeTester::Lacunary<BScal, BVect> & lat);

        /**
         * \name Sets and gets the values of <tt>m_rho</tt> and <tt>m_lossRho</tt>.
         *
         * @{
         */
        Int getRho() const { return m_rho; }
        Int getLossRho() const { return m_lossRho; }
        void setRho (const Int & val) { m_rho = val; }
        void setLossRho (const Int & val) { m_lossRho = val; }
        /*
         * @}
         */

        /**
         * Returns a non-mutable copy of the multipliers (coefficients) of the
         * MRG.
         */
        const MVect & getCoef() const { return m_aCoef; }

        /**
         * Returns the vector of multipliers \f$A\f$ as a string.
         */
        std::string toStringCoef() const;
        /**
         * The components of the lattice when it is built out of more than one
         * component. When there is only one component, it is unused as the
         * parameters are the same as above.
         */
        std::vector<MRGComponent<Int> *> comp;

      protected:

        /**
         * Initializes a square matrix of order \f$k\f$. This initial matrix 
         * contains a system of generators for the given group of states.
         */
        void initStates ();

        /**
         * Initializes some of the local variables.
         */
        void init();

        /**
         * Initializes this object when the lattice type is `ORBIT`.
         */
        void initOrbit();
        void insertion (BVect & Sta);
        void lemme2 (BVect & Sta);

        /**
         * For debugging purposes.
         */
        void trace (char* msg, int d);

        /**
         * Increments the basis by 1 in case of non-lacunary indices.
         */
        virtual void incDimBasis ();

        /**
         * Increments the basis by 1 in case of lacunary indices.
         * Uses the method described in the article: P. L'Ecuyer and R. Couture, 
         * "An Implementation of the Lattice and Spectral Tests for Multiple 
         * Recursive Linear Random Number Generators", INFORMS Journal on 
         * Computing, 9, 2 (1997), page 206--217. Section 3, "Lacunary indices".
         */
        void incDimLaBasis (int);

        /**
         * Builds the basis of the MRG recurrence in case of non-lacunary
         * indices.
         */
        void buildNaBasis (int d);

        /**
         * Builds the basis of the MRG recurrence in case of lacunary indices.
         */
        void buildLaBasis (int d);

        /**
         * \name Used for the calculation of a combined MRG.
         *
         * @{
         */
        Int m_lossRho;
        Int m_rho;
        /**
         * @}
         */

        /**
         * The coefficients of the recurrence.
         */
        MVect m_aCoef;

        /**
         * Indicates which lattice or sublattice is analyzed.
         */
        LatticeType m_latType;

        /**
         * Is `true` in the case of lacunary indices, `false` otherwise.
         */
        bool m_lacunaryFlag;

        /**
         * Contains the lacunary indices when `LacunaryFlag` is `true`,
         * otherwise is undefined.
         */
        LatticeTester::Lacunary<BScal, BVect> m_lac;


        /**
         * Work variables.
         *
         * @{
         */
        Int m_t4, m_t5, m_t6, m_t7, m_t8, m_e;
        MVect m_xi;
        /**
         * @}
         */

        /**
         * \f$\clubsuit\f$ Seems to be use as working variables. To be 
         * completed. Erwan
         */
        BMat m_sta;

        /**
         * When the flag <tt>m_ip[i]</tt> is `true`, the \f$i\f$-th diagonal
         * element of matrix <tt>m_sta</tt> is non-zero (modulo \f$m\f$) and
         * divides \f$m\f$. Otherwise (when <tt>m_ip[i]</tt> is
         * <tt>false</tt>), the \f$i\f$-th line of matrix <tt>m_sta</tt> is
         * identically 0.
         */
        bool *m_ip;
    }; // End class declaration

  //===========================================================================

  /* Max order for lacunary case in this class; takes too much memory.
     For order > ORDERMAX, use subclass MRGLatticeLac instead */
#define ORDERMAX 100

  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::trace (char *mess, int d)
    {
      std::cout << "---------------------------------------------------------------"
        << "----" << std::endl;
      std::cout << mess << std::endl;
      this->setNegativeNorm();
      this->setDualNegativeNorm();
      this->updateVecNorm ();
      this->updateDualVecNorm ();
      this->write();
      //m_w.write();
      /*
         for (int i = 0; i <= d; i++)
         std::cout << " VSI " << i << "    " << m_vSI[i] << std::endl;
         std::cout << std::endl;
         for (int i = 0; i <= d; i++)
         std::cout << " WSI " << i << "    " << m_wSI[i] << std::endl;
         */
      //checkDuality ();
      d = -1;  // compiler warning
    }


  //===========================================================================

  template<typename Int>
    MRGLattice<Int>::MRGLattice(const MRGLattice<Int> &lat):
      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice (lat.m_modulo, lat.m_order,
          lat.getDim(), lat.getNorm ()), m_lac(lat.m_lac)
  {
    m_lossRho = lat.m_lossRho;
    m_rho = lat.m_rho;
    m_latType = lat.m_latType;
    m_lacunaryFlag = lat.m_lacunaryFlag;

    m_ip = new bool[this->m_order];
    m_xi.SetLength (this->m_order);
    m_aCoef.SetLength (this->m_order);
    m_sta.SetDims (this->m_order, this->m_order);

    int dim = this->getDim();
    int rmax = std::max(this->m_order, dim);
    this->m_wSI.SetDims (rmax, dim);

    int i;
    for (i = 0; i < this->m_order; i++)
      m_aCoef[i] = lat.m_aCoef[i];
    /*
       for (i = 0; i <= m_order; i++)
       m_xi[i] = lat.m_xi[i];
       for (i = 0; i <= m_order; i++)
       m_ip[i] = lat.m_ip[i];

       int j;
       for (i = 0; i <= m_order; i++)
       for (j = 0; j <= m_order; j++)
       m_sta[i][j] = lat.m_sta[i][j];

       for (i = 0; i <= maxDim; i++)
       for (j = 0; j <= maxDim; j++)
       m_wSI[i][j] = lat.m_wSI[i][j];
       */
  }


  //===========================================================================

  template<typename Int>
    MRGLattice<Int> & MRGLattice<Int>::operator= (const MRGLattice<Int> & lat)
    {
      if (this == &lat)
        return *this;
      this->m_dim = lat.m_dim;
      copyBasis(lat);
      this->m_order = lat.m_order;
      m_ip = lat.m_ip;
      //m_shift = lat.m_shift;
      return *this;
      //MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
      //copy (lat);
      //return *this;
    }


  //===========================================================================

  template<typename Int>
    MRGLattice<Int>::MRGLattice(const Int & m, const MVect & a, int maxDim, 
        int k, LatticeType lat, LatticeTester::NormType norm):
      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice(m, k, maxDim, norm)
  {
    m_latType = lat;
    m_lacunaryFlag = false;
    m_ip = new bool[1];
    init();


    for (int i = 0; i < this->m_order; i++)
      m_aCoef[i] = a[i];
  }


  //===========================================================================

  template<typename Int>
    MRGLattice<Int>::MRGLattice(const Int & m, const MVect & a, int maxDim, 
        int k, BVect & I, LatticeType lat, LatticeTester::NormType norm):
      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice (m, k, maxDim, norm), m_lac(I, maxDim), m_ip(0)
  {
    m_latType = lat;
    m_lacunaryFlag = true;
    init();

    for (int i = 0; i < this->m_order; i++)
      m_aCoef[i] = a[i];
  }

  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::init()
    {
      kill();
      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::init();
      m_xi.SetLength(this->m_order);
      m_aCoef.SetLength(this->m_order);
      if (this->m_order > ORDERMAX) {
        m_ip = new bool[1];
        m_sta.SetDims(1, 1);
      } else {
        m_ip = new bool[this->m_order];
        m_sta.SetDims(this->m_order, this->m_order);
      }
      int rmax = std::max(this->m_order, this->getDim());
      this->m_wSI.SetDims(rmax, this->getDim());

      if (m_latType == ORBIT)
        initOrbit();
    }

  //===========================================================================

  template<typename Int>
    MRGLattice<Int>::~MRGLattice ()
    {
      kill();
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::kill()
    {
      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::kill();
      if (0 != m_ip)
        delete[] m_ip;
      m_ip = 0;
      m_xi.kill();
      m_aCoef.kill();
      m_sta.kill();
      this->m_wSI.kill();
    }

  //===========================================================================

  template<typename Int>
    BScal & MRGLattice<Int>::getLac (int j)
    {
      if (isLacunary() && j <= m_lac.getSize() && j > 0)
        return m_lac.getLac(j);
      throw std::out_of_range("MRGLattice::getLac");
    }

  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::setLac(const LatticeTester::Lacunary<BScal, BVect> & lac)
    {
      m_lac = lac;
      m_lacunaryFlag = true;
    }

  //===========================================================================

  template<typename Int>
    std::string MRGLattice<Int>::toStringCoef () const
    {
      std::ostringstream out;
      out << "[ ";
      for (int i = 0; i < this->m_order; i++)
        out << m_aCoef[i] << "  ";
      out << "]";
      return out.str ();
    }

  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::buildBasis (int d)
    {
      if (m_lacunaryFlag) {
        buildLaBasis(d);
      } else {
        buildNaBasis(d);
      }

    }

  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::buildNaBasis (int d)
    // La base est construite en dimension d.
    {
      //trace( "=====================================AVANT buildNaBasis", -10);
      initStates();


      int dk = d;
      if (dk > this->m_order)
        dk = this->m_order;

      int i, j;
      for (i = 0; i < dk; i++) {
        if (m_ip[i]) {
          for (j = 0; j < dk; j++)
            this->m_basis[i][j] = m_sta[i][j];

        } else {
          for (j = 0; j < dk; j++) {
            if (i != j)
              this->m_basis[i][j] = 0;
            else
              this->m_basis[i][j] = this->m_modulo;
          }
        }
      }

      LatticeTester::CalcDual<BMat>(this->m_basis, this->m_dualbasis, dk, this->m_modulo);
      this->setDim(dk);
      if (d > this->m_order) {
        for (i = this->m_order; i < d; i++)
          incDimBasis ();
      }
      // trace( "=================================APRES buildNaBasis", -10);
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::incDim()
    {
      if (m_lacunaryFlag) {
        incDimLaBasis (this->getDim());
      } else {
        incDimBasis ();
      }
      //   write (1);
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::incDimBasis()
    // x_n = (a_1 x_{n-1} + a_2 x_{n-2} +...+ a_k x_{n-k}) mod m. On a Dim >= Order.
    {
      // trace( "=================================AVANT incDimBasis", -10);

      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::incDim();


      const int dim = this->getDim();
      this->m_vSI.resize(dim, dim);
      this->m_wSI.resize(dim, dim);
      //m_basis.setDim(dim);
      //m_w.setDim(dim);
      //write();

      for (int i = 0; i < dim; i++) {
        clear (this->m_vSI[0][i]);
        for (int j = 0; j < this->m_order; j++) {
          conv (this->m_t1, this->m_basis[i][dim - j - 2]);
          this->m_t1 = this->m_t1 * m_aCoef[j];
          this->m_vSI[0][i] = this->m_vSI[0][i] + this->m_t1;
        }
        LatticeTester::Modulo (this->m_vSI[0][i], this->m_modulo, this->m_vSI[0][i]);
        this->m_basis[i][dim-1] = this->m_vSI[0][i];
      }

      for (int i = 0; i < dim; i++)
        this->m_basis[dim-1][i] = 0;
      this->m_basis[dim-1][dim-1] = this->m_modulo;

      for (int i = 0; i < dim-1; i++)
        this->m_dualbasis[i][dim-1] = 0;
      this->m_dualbasis[dim-1][dim-1] = 1;

      for (int j = 0; j < dim-1; j++) {

        clear (this->m_t1);
        for (int i = 0; i < dim-1; i++) {
          this->m_t2 = this->m_dualbasis[i][j];
          this->m_t2 *= this->m_vSI[0][i];
          this->m_t1 -= this->m_t2;
        }
        LatticeTester::Quotient (this->m_t1, this->m_modulo, this->m_t1);
        this->m_dualbasis[dim-1][j] = this->m_t1;
      }

      this->setNegativeNorm();
      this->setDualNegativeNorm();
      /*
         if (!checkDuality ())
         MyExit (1, "BUG");
         */
      // trace("=================================APRES incDimBasis", -10);
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::buildLaBasis (int d)
    {

      // NOT USED, see: MRGLatticeLac::buildBasis

      if (this->m_order > ORDERMAX)
        LatticeTester::MyExit (1, "MRGLattice::buildLaBasis:   k > ORDERMAX");

      initStates();
      int IMax = m_lac.getSize();

      MVect b;
      b.SetLength(this->m_order+1);
      LatticeTester::Invert(m_aCoef, b, this->m_order);

      // b is the characteristic polynomial
      PolyPE::setM (this->m_modulo);
      PolyPE::setF(b);
      PolyPE pol;
      int ord = 0;

      // Construction d'un systeme generateur modulo m.
      for (int k = 0; k < IMax; k++) {
        // pour chaque indice lacunaire
        conv (m_e, m_lac[k]);

        // x^m_e Mod f(x) Mod m
        pol.powerMod(m_e);
        pol.toVector (m_xi);

        ord = 0;
        for (int i = 1; i <= this->m_order; i++) {
          if (m_ip[i]) {
            ++ord;
            m_t5 = 0;
            for (int j = 1; j <= this->m_order; j++)
              m_t5 += m_sta[i][j] * m_xi[j - 1];
            this->m_wSI[ord][k] = m_t5;
          }
        }
      }

      /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire 
       * les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur 
       * V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes 
       * avec m sur la diagonale). 
       * */
      LatticeTester::Triangularization<BMat>(this->m_wSI, this->m_vSI, ord, IMax, 
          this->m_modulo);
      LatticeTester::CalcDual<BMat>(this->m_vSI, this->m_wSI, IMax, this->m_modulo);

      // Construire la base de dimension 1
      this->m_basis[0][0] = this->m_vSI[0][0];
      this->m_dualbasis[0][0] = this->m_wSI[0][0];
      this->setDim(1);

      this->setNegativeNorm();
      this->setDualNegativeNorm();

      for (int i = 2; i <= d; i++)
        incDimLaBasis (IMax);

      // for debugging
      // trace("ESPION_1", 1);
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::incDimLaBasis(int IMax)
    {

      LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal>::incDim();
      const int dim = this->getDim (); // new dimension (dim++)

      /*
         if (dim >= IMax) {
         MyExit (0,
         "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
         }
         */

      BVect tempLineBasis (dim);
      BVect tempColBasis (dim);

      for (int i = 0; i < dim-1; i++) {

        // tempLineBasis <- m_basis[i]
        for (int k = 0; k < dim-1; k++)
          tempLineBasis[k] = this->m_basis[i][k];

        for (int i1 = 0; i1 < dim-1; i1++) {

          BScal tempScalDual;
          LatticeTester::ProdScal<Int> (tempLineBasis, this->m_wSI[i1], dim, 
              tempScalDual);
          LatticeTester::Quotient (tempScalDual, this->m_modulo, tempScalDual);
          this->m_t1 = tempScalDual * this->m_vSI[i1][dim - 1];
          tempColBasis[i] += this->m_t1;
        }
        LatticeTester::Modulo (tempColBasis[i], this->m_modulo, tempColBasis[i]);
        this->m_basis[i][dim-1] = tempColBasis[i];
      }

      for (int j = 0; j < dim-1; j++)
        this->m_basis[dim - 1][j] = 0;
      this->m_basis[dim -1][dim - 1] = this->m_vSI[dim -1][dim - 1];

      for (int i = 0; i < dim-1; i++)
        this->m_dualbasis[i][dim - 1] = 0;

      for (int j = 0; j < dim-1; j++) {

        BScal tempScalDualBis;

        for (int i = 0; i < dim-1; i++) {
          this->m_t1 = this->m_dualbasis[i][j];
          this->m_t1 *= tempColBasis[i];
          tempScalDualBis += this->m_t1;
        }
        if (tempScalDualBis != 0)
          tempScalDualBis = -tempScalDualBis;

        LatticeTester::Quotient (tempScalDualBis, this->m_vSI[dim - 1][dim - 1], 
            tempScalDualBis);
        this->m_dualbasis[dim - 1][j] = tempScalDualBis;
      }

      LatticeTester::Quotient (this->m_modulo, this->m_vSI[dim - 1][dim - 1], this->m_t1);
      this->m_dualbasis[dim - 1][dim - 1] = this->m_t1;

      this->setNegativeNorm ();
      this->setDualNegativeNorm ();

    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::initStates ()
    /*
     * Initialise la matrice carrée Sta. La matrice Sta est d'ordre égal à
     * l'ordre du générateur. Elle contient un système de générateurs pour le
     * groupe d'états considérés.
     */
    {
      BVect statmp;
      statmp.resize(this->m_order); // Stocks variables
      int maxDim = this->getDim();
      //clear (m_t2);

      if (m_latType == RECURRENT) {
        // check if a_k is relatively prime to m --> m_t1 = 1
        this->m_t1 = GCD (m_aCoef[this->m_order], this->m_modulo);
        this->m_t1 = abs(this->m_t1);
        LatticeTester::set9 (this->m_t2);
      }

      if (m_latType == FULL || m_latType == PRIMEPOWER || (this->m_t1 == this->m_t2)) {
        // m_sta is set to identity matrix
        for (int i = 0; i < this->m_order; i++) {
          for (int j = 0; j < this->m_order; j++) {
            if (i != j)
              clear (m_sta[i][j]);
            else
              LatticeTester::set9 (m_sta[i][j]);
          }
          m_ip[i] = true;
        }
        double temp;
        conv(temp, this->m_modulo);
        double lgm2 = 2.0 * LatticeTester::Lg (temp);
        this->calcLgVolDual2 (lgm2);

      } else {
        if (m_latType == ORBIT) {
          LatticeTester::MyExit (1, "case ORBIT is not finished");

          MVect InSta;
          Int inStatmp;
          InSta.SetLength (this->m_order);
          clear (statmp[this->m_order-1]);
          for (int i = 0; i < this->m_order; i++) {
            InSta[0] = m_aCoef[i] * InSta[this->m_order - i - 1];
            statmp[this->m_order-1] += InSta[0];
          }

          statmp[this->m_order-1] -= InSta[this->m_order];
          for (int i = 0; i < this->m_order; i++)
            statmp[i] = InSta[i+1] - InSta[i];
          InSta.kill();

        } else if (m_latType == RECURRENT) {
          PolyPE::setM (this->m_modulo);
          /* Je crois que la version sunos devait fonctionner correctement.
           * Je crois qu'Ajmal a créé des bugs dans la version mcs, qui se sont 
           * propagés à mcs2, xds98, ..., c++. Je n'ai pas réussi à trouver 
           * l'erreur: les résultats de plusieurs exemples dans sunos ne 
           * concordent pas avec les résultats des versions subséquentes. Voir 
           * l'exemple 4 dans l'article
           *    AUTHOR="P. L'Ecuyer and R. Couture",
           *    TITLE="An Implementation of the Lattice and Spectral Tests for
           *           Multiple Recursive Linear Random Number Generators"
           * A CORRIGER: comparer avec /u/lecuyer/stochas/latmrg/sunos/
           * Voir la déf du m effectif comparé au vrai m dans LATIO. 
           * Je soupçonne que cela pourrait être l'origine des erreurs.
           */

          LatticeTester::MyExit (1, "case RECURRENT ne fonctionne pas");
          printf("ESPION_RECURRENT\n");
          MVect b;
          b.SetLength(this->m_order + 1);
          LatticeTester::CopyVect (b, m_aCoef, this->m_order);
          PolyPE::reverse (b, this->m_order, 2);
          // b is the characteristic polynomial
          PolyPE::setF(b);
          PolyPE pol;

          // Must have 2^m_e > m^k to be sure to reach a recurrent state
          m_e = 3 + (int) (this->m_order * 0.5 * this->m_lgm2);
          pol.powerMod(m_e);
          pol.toVector (m_xi);

          statmp[0] = m_xi[this->m_order - 1];
          for (int i = 2; i <= this->m_order; i++) {
            // Multiplier m_xi par X et reduire mod X^k - a1 X^{k-1} - ....

            m_xi[this->m_order] = m_xi[this->m_order - 1];
            for (int j = 1; j < this->m_order; j++) {
              // Coeff. de X^{m_order-j}.
              m_xi[this->m_order - j] = m_xi[this->m_order - j - 1];
              // ********* ATTENTION: MulMod (a, b, c, d) est très différent
              // pour les types long et ZZ (voir ZZ.txt). C'est pourquoi on
              // utilise MulMod (a, b, c) même s'il est plus lent. *********
              m_xi[this->m_order - j - 1] = MulMod (m_xi[this->m_order], m_aCoef[j], 
                  this->m_modulo);
              m_xi[this->m_order - j] += m_xi[this->m_order - j - 1];
            }
            // Coeff. constant.
            m_xi[0] = MulMod (m_xi[this->m_order], m_aCoef[this->m_order], this->m_modulo);
            statmp[i] = m_xi[this->m_order - 1];
          }
        }

        for (int i = 0; i < this->m_order; i++) {
          for (int j = 0; j < this->m_order; j++)
            clear (m_sta[i][j]);
          m_ip[i] = false;
          m_sta[i][0] = statmp[i];
        }
        insertion (statmp);

        for (int k = 1; k < this->m_order; k++) {
          // On passe a l'etat suivant.
          for (int j = 0; j < this->m_order-1; j++)
            statmp[j] = m_sta[j + 1][0];
          clear (statmp[this->m_order-1]);
          for (int i = 0; i < this->m_order; i++) {
            this->m_t1 = m_aCoef[i] * m_sta[this->m_order - i + 1][0];
            statmp[this->m_order-1] += this->m_t1;
          }
          LatticeTester::Modulo(statmp[this->m_order-1], this->m_modulo, 
              statmp[this->m_order-1]);
          // On memorise l'etat suivant.
          for (int i = 0; i < this->m_order-1; i++)
            swap (m_sta[i][0], m_sta[i + 1][0]);
          m_sta[this->m_order-1][0] = statmp[this->m_order-1];
          insertion (statmp);
        }

        lemme2 (statmp);

        // Calcul de lgVolDual2
        double x;
        if (m_ip[1]) {
          conv(x, this->m_modulo / m_sta[0][0]);
          this->m_lgVolDual2[0] = 2.0 * LatticeTester::Lg (x);
        } else
          this->m_lgVolDual2[0] = 0.0;

        int rmax = std::min(this->m_order, maxDim);
        for (int r = 2; r <= rmax; r++) {
          if (m_ip[r]) {
            conv(x, this->m_modulo / m_sta[r][r]);
            this->m_lgVolDual2[r] = this->m_lgVolDual2[r-1] + 2.0 * LatticeTester::Lg (x);
          } else
            this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
        }

        for (int r = this->m_order + 1; r <= maxDim; r++)
          this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
      }
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::insertion (BVect & statmp)
    /*
     * Cette procedure insere le vecteur Sta[0] dans la matrice triangulaire
     * Sta. Si IP[i] = TRUE, l'entree diagonale sur la i-ieme ligne de Sta est
     * non-nulle (modulo m) et divise m. Sinon, la i-ieme ligne est
     * identiquement nulle. L'insertion doit preserver ces proprietes.
     * Le vecteur Sta[0] est altere au cours de l'operation.
     */
    {
      for (int j = 0; j < this->m_order; j++) {
        LatticeTester::Modulo (statmp[j], this->m_modulo, statmp[j]);
        if (!IsZero (statmp[j])) {
          if (!m_ip[j]) {
            LatticeTester::Euclide (statmp[j], this->m_modulo, this->m_t1, this->m_t2, this->m_t3, 
                m_t4, m_sta[j][j]);
            for (int i = j + 1; i < this->m_order; i++) {
              m_sta[j][i] = this->m_t1 * statmp[i];
              LatticeTester::Modulo (m_sta[j][i], this->m_modulo, m_sta[j][i]);
            }
            m_ip[j] = true;
            return;

          } else {
            LatticeTester::Euclide (m_sta[j][j], statmp[j], this->m_t1, this->m_t2, this->m_t3, 
                m_t4, m_sta[j][j]);
            clear (statmp[j]);
            for (int i = j + 1; i < this->m_order; i++) {
              m_t5 = this->m_t1 * m_sta[j][i];
              m_t6 = this->m_t2 * statmp[i];
              m_t7 = this->m_t3 * m_sta[j][i];
              m_t8 = m_t4 * statmp[i];
              m_sta[j][i] = m_t5 + m_t6;
              LatticeTester::Modulo (m_sta[j][i], this->m_modulo, m_sta[j][i]);
              statmp[i] = m_t7 + m_t8;
            }
          }
        }
      }
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::lemme2 (BVect & statmp)
    /*
     * Cette procedure suppose que la matrice Sta est triangulaire. Si
     * IP[i] = TRUE, l'entree diagonale sur la i-ieme ligne de Sta est
     * non-nulle (modulo m) et divise m. Sinon, la i-ieme ligne est
     * identiquement nulle.
     */
    {
      for (int i = 0; i < this->m_order; i++) {
        if (m_ip[i]) {
          LatticeTester::Quotient (this->m_modulo, m_sta[i][i], this->m_t1);
          this->m_t1 = abs (this->m_t1);
          if (this->m_t1 < this->m_modulo) {
            for (int j = 0; j < i; j++)
              statmp[j] = m_sta[i][j];
            clear (m_sta[0][i]);
            for (int j = i + 1; j < this->m_order; j++)
              statmp[j] = this->m_t1 * m_sta[i][j];
            insertion (statmp);
          }
        }
      }
    }


  //===========================================================================

  template<typename Int>
    void MRGLattice<Int>::initOrbit()
    {
      LatticeTester::MyExit (1, "MRGLattice::initOrbit n\'est pas terminée.");

      /*
         for (int j = 0; j < J; j++) {
         for (int i = 1; i <= k; i++) {
         if (j == 0)
         clear(InSta[i]);
         for (int i1 = 1; i1 <= kj; i1++) {
         Multiply (aj [j,i1], VectSup [i-i1], SupT3);
         Add (SupT3, VectSup [i], VectSup [i])
         }
         Multiply (nj [j], VectSup [i], SupT3);
         Add (InSta [i], SupT3, InSta [i]);
         if (j == J - 1) {  Modulo (InSta [i], mm, InSta [i])  }
         }
         }
         */
    }


} // End namespace LatMRG
#endif

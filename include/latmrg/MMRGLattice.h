#ifndef MMRGLATTICE_H
#define MMRGLATTICE_H
//#include "latticetester/Types.h"
//#include "latticetester/Const.h"
//#include "latticetester/Lacunary.h"
//#include "latticetester/IntLattice.h"

//#include "latmrg/Const.h"

//#include <string>


namespace LatMRG {

  /**
   * This class implements lattice basis built from M-MRG (matrix multiple
   * recursive linear congruential generators). One must first call the 
   * constructor with a given congruence modulus \f$m\f$, a given generator
   * matrix for the recurrence, and a maximal dimension for the basis. One must 
   * then build the lattice basis associated to the generator matrix for a
   * given dimension. Each MMRG is defined by a generator matrix \f$A\f$. This 
   * MMRG satisfies the recurrence
   * \f[
   *   X_n = A X_{n-1} \mod m.
   * \f]
   */
  template<typename Int>
    class MMRGLattice: public LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> {
      public:

        /**
         * Constructor with modulus of congruence \f$m\f$, generator matrix 
         * \f$A\f$, dimension of generator matrix \f$r\f$, maximal dimension 
         * `MaxDim`, and lattice type `Latt`. Vectors and (square) matrices of 
         * the basis have maximal dimension `maxDim`, and the indices of 
         * vectors and matrices vary from dimension 0 to `maxDim`-1. The norm 
         * to be used for the basis vectors is `norm`.
         */
        MMRGLattice (const Int & m, const MMat & A, int maxDim, int r,
            LatticeTester::NormType norm = LatticeTester::L2NORM,
            LatticeType lat = FULL);

        /**
         * As in the constructor above but the basis is built for the lacunary
         * indices `lac`.
         */
        //PW_TODO à faire plus tard
        MMRGLattice (const Int & m, const MMat & A, int maxDim, int r,
            LacunaryType & lacunaryType, BVect & lac, 
            LatticeTester::NormType norm = LatticeTester::L2NORM, 
            LatticeType lat = FULL);

        /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        MMRGLattice (const MMRGLattice<Int> & Lat);

        /**
         * Destructor.
         */
        ~MMRGLattice();

        /**
         * Cleans and releases memory used by this object.
         */
        void kill();

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set equal to <tt>Lat</tt>’s current dimension.
         */
        MMRGLattice<Int> & operator= (const MMRGLattice<Int> & Lat);

        /**
         * Returns the \f$j\f$-th lacunary index.
         */
        BScal & getLac (int j);

        /**
         * Sets the lacunary indices for this lattice to `lat`.
         */
        virtual void setLac (const LatticeTester::Lacunary<BScal, BVect> & lat);

        /**
         * Returns the generator matrix \f$A\f$ as a string.
         */
        std::string toStringGeneratorMatrix() const;

        /**
         * Builds the basis in dimension \f$d\f$.
         */
        virtual void buildBasis (int d);

        /**
         * Increments the dimension of the basis by 1 by calling either
         * `incDimBasis` or `incDimLaBasis`.
         */
        void incDim();

        /**
         * Returns `true` for the case of lacunary indices, returns `false` for
         * non-lacunary indices.
         */
        bool isLacunary() const { return m_lacunaryFlag; }

        /**
         * Returns a non-mutable copy of the generator matrix of the MMRG
         */
        const MMat & getGeneratorMatrix() const { return m_A; }


        //PW_TODO ici temporairement, à déplacer dans latticetester/Util.h
        void getSubLine(MVect & vec, MMat& B, int lign, int jMin, int jMax);


      protected:

        /**
         * Initializes some of the local variables.
         */
        void init();

        /**
         * Builds the basis of the MMRG recurrence in case of non-lacunary
         * indices.
         */
        void buildNonLacunaryBasis (int dimension);

        /**
         * Builds the basis of the MMRG recurrence in case of lacunary
         * indices.
         */
        void buildLacunaryBasis (int dimension);

        /**
         * Increments the basis by 1 in case of non-lacunary indices.
         */ 
        void incrementDimNonLacunaryBasis ();

        /**
         * Increments the basis by 1 in case of lacunary indices.
         */
        void incrementDimLacunaryBasis(int Imax);

        /**
         * The generator matrix of the recurrence.
         */
        MMat m_A;

        /**
         * Indicates which lattice or sublattice is analyzed.
         */
        LatticeType m_latType;

        /**
         * Is `true` in the case of lacunary indices, `false` otherwise.
         */
        bool m_lacunaryFlag;

        /**
         * Type of the lacunary projection selected.
         */
        LacunaryType m_lacunaryType;
        /**
         * Contains the lacunary indices when `LacunaryFlag` is `true`,
         * otherwise is undefined.
         */
        LatticeTester::Lacunary<BScal, BVect> m_lac;

        /**
         * Contains the number of lacunary indices
         */ 
        int m_numberLacIndices;

        /**
         * Matrix used for lacunary indices
         */
        BMat m_B;

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
         * \f$\clubsuit\f$ Seems to be use as working variables. 
         * To be completed. Erwan
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

  /* Max order for lacunary case in this class; takes too much memory.
     For order > ORDERMAX, use subclass MRGLatticeLac instead */
  //PW_TODO à voir plus tard avec lacunary
#define ORDERMAX 100

  //===========================================================================

  template<typename Int>
    MMRGLattice<Int>::MMRGLattice(const Int & m, const MMat & A, int maxDim, 
        int r, LatticeTester::NormType norm, LatticeType lat):
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice(m, r, maxDim, norm)
  {
    m_A = A;
    m_latType = lat;
    m_lacunaryFlag = false;
    m_lacunaryType = NONE;
    //m_ip = new bool[1];
    init();
    //PW_TODO attention m_ip aussi initialisé dans init()
  }


  //===========================================================================

  template<typename Int>
    MMRGLattice<Int>::MMRGLattice(const Int & m, const MMat & A, int maxDim, 
        int r, LacunaryType & lacunaryType, BVect & lac, 
        LatticeTester::NormType norm, LatticeType lat):
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice (m, r, maxDim, norm)
      //m_lac(lac, r)
      //PW_TODO r ou maxDim?
  {
    m_ip=0;
    m_lac = LatticeTester::Lacunary<BScal, BVect>(lac, maxDim);
    m_A = A;
    m_latType = lat;
    m_lacunaryFlag = true;
    m_lacunaryType = lacunaryType;
    m_numberLacIndices = m_lac.getSize();
    init();
  }


  //===========================================================================

  template<typename Int>
    MMRGLattice<Int>::MMRGLattice(const MMRGLattice & lat):
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::IntLattice (lat.m_modulo, lat.m_order,
          lat.getDim(), lat.getNorm ()), m_lac(lat.m_lac)
  {
    m_A = lat.m_A;
    m_latType = lat.m_latType;
    m_lacunaryFlag = lat.m_lacunaryFlag;

    m_ip = new bool[m_order];
    m_xi.SetLength (m_order);
    m_A.SetDims (m_order, m_order);
    m_sta.SetDims (m_order, m_order);

    int dim = getDim();
    int rmax = std::max(m_order, dim);
    m_wSI.SetDims (rmax, dim);

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
    //PW_TODO : ça doit vraiment rester commenté ?
  }

  //===========================================================================

  template<typename Int>
    MMRGLattice<Int>::~MMRGLattice ()
    {
      kill();
    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::kill()
    {
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::kill();
      if (0 != m_ip)
        delete[] m_ip;
      m_ip = 0;
      m_xi.kill();

      // PW_TODO : methode kill fonctionne sur une matrice ?
      //m_aCoef.kill();
      //m_A.kill();

      // PW_TODO à quoi ça sert de killer ça ?
      m_sta.kill();
      m_wSI.kill();
      //m_vSI.kill(); ??
    }

  //===========================================================================

  template<typename Int>
    MMRGLattice<Int> & MMRGLattice<Int>::operator= 
    (const MMRGLattice<Int> & lat)
    {
      if (this == &lat)
        return *this;
      m_dim = lat.m_dim;
      copyBasis(lat);
      m_order = lat.m_order;
      m_ip = lat.m_ip;
      //m_shift = lat.m_shift;
      return *this;
      //MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
      //copy (lat);
      //return *this;
    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::init()
    {
      kill(); //PW_TODO : wzf ?
      IntLattice::init();
      m_xi.SetLength(m_order);
      m_A.SetDims(m_order, m_order);
      if (m_order > ORDERMAX) {
        m_ip = new bool[1];
        m_sta.SetDims(1, 1);
      } else {
        m_ip = new bool[m_order];
        m_sta.SetDims(m_order, m_order);
      }
      int rmax = std::max(m_order, getDim());
      m_wSI.SetDims(rmax, getDim());

      double temp;
      NTL::conv(temp, m_modulo);
      double lgm2 = 2.0 * LatticeTester::Lg (temp);
      calcLgVolDual2 (lgm2);
      //if (m_latType == ORBIT)
      //   initOrbit();
      //PW_TODO
    }

  //===========================================================================

  template<typename Int>
    BScal & MMRGLattice<Int>::getLac (int j)
    {
      if (isLacunary() && j <= m_lac.getSize() && j > 0)
        return m_lac.getLac(j);
      throw std::out_of_range("MMRGLattice::getLac");
    }


  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::setLac(const LatticeTester::Lacunary<BScal, BVect> & lac)
    {
      m_lac = lac;
      m_lacunaryFlag = true;
    }


  //===========================================================================

  //PW_TODO merdasse à tester

  template<typename Int>
    std::string MMRGLattice<Int>::toStringGeneratorMatrix () const
    {
      std::ostringstream out;
      out << "[";
      for (int i = 0; i < m_order; i++) {
        out << "[";
        for (int j = 0; j < (m_order-1); j++) {
          out << m_A[i][j] << " ";
        }
        if (i == (m_order-1))
          out << m_A[i][m_order-1] << "]";
        else
          out << m_A[i][m_order-1] << "]" << std::endl;
      }
      out << "]" << std::endl;

      return out.str ();
    }


  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::buildBasis (int d)
    {
      if (m_lacunaryFlag)
        buildLacunaryBasis(d);
      else
        buildNonLacunaryBasis(d);
    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::buildNonLacunaryBasis (int dimension)
    // a basis is built in dimension d

    {
      setDim(dimension);
      int sizeA = getOrder();
      m_basis.resize(dimension, dimension);
      m_vecNorm.resize(dimension);
      m_dualvecNorm.resize(dimension);
      setNegativeNorm ();
      setDualNegativeNorm ();

      // filling in the diagonal of m_basis
      for (int i = 0; i < sizeA; i++)
        m_basis[i][i] = 1;
      for (int i = sizeA; i < dimension; i++)
        m_basis[i][i] = m_modulo;

      // using genrator matrix A to complete the first lines of m_basis
      // with values generated by the recurrence
      NTL::ZZ_p::init(m_modulo);
      MMatP temp;
      temp.SetDims(sizeA, sizeA);
      for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

      int maxIter = floor(dimension/sizeA);

      for (int k = 1; k < maxIter+1; k++) {
        // calculation of transpose(A^k)
        temp *= NTL::conv<MMatP>(transpose(m_A));

        if (k == maxIter) { // we completed the end of m_basis matrix
          int residu = dimension - maxIter * sizeA;
          for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < residu; j ++)
              m_basis[i][k*sizeA +j] = NTL::conv<Int>(temp[i][j]);
          }
        } else {
          for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < sizeA; j ++)
              m_basis[i][k*sizeA +j] = NTL::conv<Int>(temp[i][j]);
          }
        }
      }

      // we create the dual lattice associated
      m_dualbasis.resize(dimension, dimension);
      LatticeTester::CalcDual<BMat>(m_basis, m_dualbasis, dimension, m_modulo);

      if (!checkDuality())
        LatticeTester::MyExit (1, "BUG in MMRGLattice::buildNonLacunaryBasis");
    }


  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::getSubLine(MVect & vec, MMat& B, int lign, int jMin,
        int jMax)
    {
      // both jMin and jMax are included
      vec.resize(jMax-jMin+1);
      for (int i = 0; i < (jMax-jMin+1); i++)
        vec[i] = B[lign][jMin+i];
    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::buildLacunaryBasis (int dimension)
    {
      int sizeA = getOrder();
      m_vecNorm.resize(dimension);
      m_dualvecNorm.resize(dimension);

      int maxIndiceLac = NTL::conv<int>(m_lac[m_lac.getSize()-1]); 

      // building the complete basis until: dimension = max lacunary indice
      //-----------------------------------------------------------------------
      BMat tempBasis;
      //PW_TODO
      tempBasis.resize(m_order, maxIndiceLac+1);
      // PW_TODO meilleurs size à trouver
      // +1 because lacunary indices start at 0

      // filling in the diagonal of tempBasis
      for (int i = 0; i < sizeA; i++)
        tempBasis[i][i] = 1;
      //for (int i = sizeA; i < maxIndiceLac; i++)
      //   tempBasis[i][i] = m_modulo;

      // using genrator matrix A to complete the first lines of tempBasis
      // with values generated by the recurrence
      NTL::ZZ_p::init(m_modulo);
      MMatP temp;
      temp.SetDims(sizeA, sizeA);
      for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

      int maxIter = floor((maxIndiceLac+1)/sizeA);

      for (int k = 1; k < maxIter+1; k++) {
        // calculation of transpose(A^k)
        temp *= NTL::conv<MMatP>(transpose(m_A));

        if (k == maxIter) { // we completed the end of m_basis matrix
          int residu = (maxIndiceLac+1) - maxIter * sizeA;
          for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < residu; j ++)
              tempBasis[i][k*sizeA +j] = NTL::conv<Int>(temp[i][j]);
          }
        } else {
          for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < sizeA; j ++)
              tempBasis[i][k*sizeA +j] = NTL::conv<Int>(temp[i][j]);
          }
        }
      }

      // projecting over the columns of interest (lacunary indices)
      //-----------------------------------------------------------------------
      m_wSI.resize(std::max(m_order, m_numberLacIndices), m_numberLacIndices);
      m_vSI.resize(std::max(m_order, m_numberLacIndices), m_numberLacIndices);
      // PW_TODO meilleurs size à trouver

      for (int j = 0; j < m_numberLacIndices; j++) {
        for (int i = 0; i < m_order; i++)
          m_wSI[i][j] = tempBasis[ i ][ NTL::conv<int>(m_lac[j]) ];
      }

      /*
         std::cout << "tempBasis = \n" << tempBasis << std::endl;
         std::cout << "projection = \n" << m_wSI << std::endl;
         std::cout << "\n******************************************\n" << std::endl;
         std::cout << "\nAVANT TRIANGULARIZATION" << std::endl;
         std::cout << "m_vSI = \n" << m_vSI << std::endl;
         std::cout << "m_wSI = \n" << m_wSI << std::endl;
         */

      // transforming this generating familly into a basis of the lattice
      //-----------------------------------------------------------------------
      LatticeTester::Triangularization <BMat> (m_wSI, m_vSI, m_order,
          m_numberLacIndices, m_modulo);

      //std::cout << "\nAPRES TRIANGULARIZATION, AVANT CALCDUAL" << std::endl;
      //std::cout << "m_vSI = \n" << m_vSI << std::endl;
      //std::cout << "m_wSI = \n" << m_wSI << std::endl;

      LatticeTester::CalcDual <BMat> (m_vSI, m_wSI, m_numberLacIndices,
          m_modulo);

      /*
         std::cout << "\nAPRES CALCDUAL" << std::endl;
         std::cout << "m_vSI = \n" << m_vSI << std::endl;
         std::cout << "m_wSI = \n" << m_wSI << std::endl;
         std::cout << "\n******************************************\n" << std::endl;
         */

      //building the basis in dimension 1
      m_basis[0][0] = m_vSI[0][0];
      m_dualbasis[0][0] = m_wSI[0][0];
      setDim (1);

      setNegativeNorm();
      setDualNegativeNorm();


      for (int i = 1; i < dimension; i++)
        incDim ();

      if (!checkDuality())
        LatticeTester::MyExit (1, "BUG in MMRGLattice::buildNonLacunaryBasis");

    }


  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::incDim()
    {
      if (m_lacunaryFlag)
        incrementDimLacunaryBasis(m_numberLacIndices);
      else
        incrementDimNonLacunaryBasis();
    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::incrementDimNonLacunaryBasis()
    // X_n = A X_{n-1} mod m. We have: dimension >= order.
    {
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::incDim();
      int newDimension = getDim();
      int sizeA = getOrder();

      // ************* update of the primal lattice *************
      //  - we add a new coordinate to each vector v_i, this value being 
      //  determined by the MMRG recurrence (even if the original vectors have 
      //  been transformed linearly and we must apply the same transformations 
      //  to their last coordinates).
      //  - we add an extra vector (0,..., 0, m) to complete this dimension
      //    increased basis.

      // we compute the number of steps required to reach the current state of
      // the generator for the considered dimension.
      int n = floor((newDimension-1) / sizeA);
      NTL::ZZ_p::init(m_modulo);
      MMatP temp;
      temp.SetDims(sizeA, sizeA);
      for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;
      for (int k = 0; k < n; k++)
        temp *= NTL::conv<MMatP>(transpose(m_A));
      // PW_TODO : could be useful to keep A^k in memory to shorten computation

      // update of the new v_i coordinates using the *temp* matrix and the 
      // first coefficients of each line (can be seen as a seed vector). So 
      // this *temp* matrix multiplied by this seed vector gives us the next 
      // values generated by the MMRG for the considered dimension.
      MVect initialState;
      for (int i = 0; i < (newDimension-1); i++) {
        getSubLine(initialState, m_basis, i, 0, sizeA-1);
        initialState = NTL::conv<MVect>( 
            transpose(temp) * NTL::conv<MVectP>(initialState));
        m_basis[i][newDimension-1] = initialState[newDimension - n*sizeA -1];
      }
      m_basis[newDimension-1][newDimension-1] = m_modulo;


      // ************* update of the dual basis *************
      //  - we add a new 0 coordinate to each vector w_i in the dual basis.
      //  - for the new last line of the matrix, we add an extra vector,
      //    as described in L'Ecuyer's paper
      // PW_TODO : citer référence "Guide LatTester"

      MVect lastLine;
      lastLine.SetLength(newDimension);

      for (int i = 0; i < (newDimension-1); i++) {
        NTL::matrix_row<const BMat> row(m_dualbasis, i);
        lastLine -= m_basis[i][newDimension-1] * row;
      }
      for (int i = 0; i < (newDimension-1); i++)
        m_dualbasis[newDimension-1][i] = lastLine[i] / m_modulo;
      m_dualbasis[newDimension-1][newDimension-1] = 1;

      setNegativeNorm();
      setDualNegativeNorm();

      if (!checkDuality())
        LatticeTester::MyExit (1, "BUG in MMRGLattice::incrementDimBasis");

    }

  //===========================================================================

  template<typename Int>
    void MMRGLattice<Int>::incrementDimLacunaryBasis(int Imax)
    {
      IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal>::incDim();
      const int dim = getDim (); // new dimension (dim++)

      //PW_TODO
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
          tempLineBasis[k] = m_basis[i][k];

        // v[i] -> VSI[0].
        // for (int j = 0; j < dim; j++)
        //    m_vSI[0][j] = m_basis[i][j];

        //clear (m_vSI[i][0]);

        for (int i1 = 0; i1 < dim-1; i1++) {

          BScal tempScalDual;

          LatticeTester::ProdScal<Int> (tempLineBasis, m_wSI[i1], dim,
              tempScalDual);
          LatticeTester::Quotient (tempScalDual, m_modulo, tempScalDual);
          m_t1 = tempScalDual * m_vSI[i1][dim - 1];
          tempColBasis[i] += m_t1;
        }
        LatticeTester::Modulo (tempColBasis[i], m_modulo, tempColBasis[i]);
        m_basis[i][dim-1] = tempColBasis[i];
      }

      for (int j = 0; j < dim-1; j++)
        m_basis[dim - 1][j] = 0;
      m_basis[dim -1][dim - 1] = m_vSI[dim -1][dim - 1];

      for (int i = 0; i < dim-1; i++)
        m_dualbasis[i][dim - 1] = 0;

      for (int j = 0; j < dim-1; j++) {

        //clear (m_wSI[0][j]);
        BScal tempScalDualBis;

        for (int i = 0; i < dim-1; i++) {
          m_t1 = m_dualbasis[i][j];
          m_t1 *= tempColBasis[i];
          tempScalDualBis += m_t1;
        }
        if (tempScalDualBis != 0)
          tempScalDualBis = -tempScalDualBis;

        LatticeTester::Quotient (tempScalDualBis, m_vSI[dim - 1][dim - 1],
            tempScalDualBis);
        m_dualbasis[dim - 1][j] = tempScalDualBis;
      }

      LatticeTester::Quotient (m_modulo, m_vSI[dim - 1][dim - 1], m_t1);
      m_dualbasis[dim - 1][dim - 1] = m_t1;

      //setDim (dim + 1);
      setNegativeNorm ();
      setDualNegativeNorm ();

    }

} // End namespace LatMRG
#endif

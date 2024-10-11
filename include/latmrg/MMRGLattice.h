#ifndef LATMRG_MMRGLATTICE_H
#define LATMRG_MMRGLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Types.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/MRGLattice.h"
#include "latmrg/FlexModInt.h"
//#include "Primitivity.h"

namespace LatMRG {

  /**
   * \todo Make this class not use/compute the dual if not asked.
   *
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
 template<typename Int, typename Real>
class MMRGLattice: public LatticeTester::IntLatticeExt<Int, Real> {
   public:
   // Parent:
   // MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

      /**
       * Constructor with modulus of congruence \f$m\f$, order of the recurrence
       * \f$k\f$, multipliers in `aa`, and maximal dimension `maxDim`.
       * The length of basis vectors is computed with `norm`.
       * The basis is built for the lacunary indices in `lac`.
       * The vector `aa` must have k+1 components with `a[j]`=\f$a_j\f$.
       */
      MMRGLattice(const Int & m, const IntMat & A, int maxDim,
            LatticeTester::NormType norm = LatticeTester::L2NORM,
            LatticeType lat = FULL);
   
      /**
       * As in the constructor above but the basis is built for the lacunary
       * indices `lac`.
      */
      //PW_TODO à faire plus tard
      MMRGLattice (const Int & m, const IntMat & A, int maxDim,
         LacunaryType & lacunaryType, IntVec & lac,
         LatticeTester::NormType norm = LatticeTester::L2NORM,
         LatticeType lat = FULL);

       /**
         * Copy constructor. The maximal dimension of the created basis is set
         * equal to <tt>Lat</tt>’s current dimension.
         */
        MMRGLattice (const MMRGLattice<Int, Real> & Lat);

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
        MMRGLattice<Int, Real> & operator= (const MMRGLattice<Int, Real> & Lat);

       /**
         * Returns the \f$j\f$-th lacunary index.
         */
        Int & getLac (int j);

        /**
         * Sets the lacunary indices for this lattice to `lat`.
         */
        virtual void setLac (const LatticeTester::Lacunary<Int> & lat);

        /**
         * Returns the generator matrix \f$A\f$ as a string.
         */
        std::string toString() const;

        /**
         * Builds the basis in dimension \f$d\f$.
         */
        void buildBasis (int64_t d) override;

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
        const IntMat & getGeneratorMatrix() const { return m_A; }


        //PW_TODO ici temporairement, à déplacer dans latticetester/Util.h
        void getSubLine(IntVec & vec, IntMat& B, int lign, int jMin, int jMax);

  

      protected:

        /**
         * Initializes some of the local variables.
         */
        void init();

        /**
         * Builds the basis of the MMRG recurrence in case of non-lacunary
         * indices.
         */
        void buildNonLacunaryBasis (int d);

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
        IntMat m_A;

        /**
         * Basis that is always of a dimension that is a multiple of the order
         * and that is always of a greater or equal dimension than the current
         * basis. This is used by `buildBasis()` and `incDim()` so that powers
         * of `A` do not have to be computed every time.
         * */
        IntMat m_basis_max;

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
        LatticeTester::Lacunary<Int> m_lac;

        /**
         * Contains the number of lacunary indices
         */
        int m_numberLacIndices;

        /**
         * Matrix used for lacunary indices
         */
        IntMat m_B;

        /**
         * Work variables.
         *
         * @{
         */
        //IntVec m_xi;
        /**
         * @}
         */

        /**
         * \f$\clubsuit\f$ Seems to be use as working variables.
         * To be completed. Erwan
         */
        //IntMat m_sta;

        /**
         * When the flag <tt>m_ip[i]</tt> is `true`, the \f$i\f$-th diagonal
         * element of matrix <tt>m_sta</tt> is non-zero (modulo \f$m\f$) and
         * divides \f$m\f$. Otherwise (when <tt>m_ip[i]</tt> is
         * <tt>false</tt>), the \f$i\f$-th line of matrix <tt>m_sta</tt> is
         * identically 0.
         */
        ///bool *m_ip;
    }; // End class declaration

  /* Max order for lacunary case in this class; takes too much memory.
     For order > ORDERMAX, use subclass MRGLatticeLac instead */
  //PW_TODO à voir plus tard avec lacunary
#define ORDERMAX 100
  //===========================================================================

  template<typename Int, typename Real>
  MMRGLattice<Int, Real>::MMRGLattice(const Int & m, const IntMat & A, int maxDim,
      LatticeTester::NormType norm, LatticeType lat) : IntLatticeExt<Int, Real>(m, maxDim, norm) {
    m_A = A;
    m_latType = lat;
    m_lacunaryFlag = false;
    m_lacunaryType = NONE;
    init(); 
   }
   
  //===========================================================================

  template<typename Int, typename Real>
    MMRGLattice<Int, Real>::MMRGLattice(const Int & m, const IntMat & A, int maxDim,
        LacunaryType & lacunaryType, IntVec & lac,
        LatticeTester::NormType norm, LatticeType lat):
      IntLatticeExt<Int, Real> (m, maxDim, norm)
      //m_lac(lac, r)
      //PW_TODO r ou maxDim?
  {
    //m_ip=0;
    m_lac = LatticeTester::Lacunary<Int>(lac, maxDim);
    m_A = A;
    m_latType = lat;
    m_lacunaryFlag = true;
    m_lacunaryType = lacunaryType;
    m_numberLacIndices = m_lac.getSize();
    init(); 
  }
  
  //===========================================================================

  template<typename Int, typename Real>
    MMRGLattice<Int, Real>::MMRGLattice(const MMRGLattice & lat):
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.m_order,
          lat.getDim(), true, lat.getNorm ()), m_lac(lat.m_lac)
  {
    m_A = lat.m_A;
    m_latType = lat.m_latType;
    m_lacunaryFlag = lat.m_lacunaryFlag;

    //m_ip = new bool[this->m_order];
    //m_xi.SetLength (this->m_order);
    m_A.SetDims (this->m_order, this->m_order);
    //m_sta.SetDims (this->m_order, this->m_order);

    int dim = this->getDim();
    int rmax = std::max(this->m_order, dim);
    this->m_wSI.SetDims (rmax, dim);

    /*
       for (i = 0; i <= this->m_order; i++)
       m_xi[i] = lat.m_xi[i];
       for (i = 0; i <= this->m_order; i++)
       m_ip[i] = lat.m_ip[i];

       int j;
       for (i = 0; i <= this->m_order; i++)
       for (j = 0; j <= this->m_order; j++)
       m_sta[i][j] = lat.m_sta[i][j];

       for (i = 0; i <= maxDim; i++)
       for (j = 0; j <= maxDim; j++)
       this->m_wSI[i][j] = lat.m_wSI[i][j];
       */
    //PW_TODO : ça doit vraiment rester commenté ?
  }

  //===========================================================================


  template<typename Int, typename Real>
    MMRGLattice<Int, Real>::~MMRGLattice ()
    {
       kill();
    }
  
  //===========================================================================
    
    template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::kill()
    {
      //if (0 != m_ip)
      //  delete[] m_ip;
      //m_ip = 0;
      //m_xi.kill();

      //m_sta.kill();
      //this->m_wSI.kill();
    }
    
  //===========================================================================

  template<typename Int, typename Real>
    MMRGLattice<Int, Real> & MMRGLattice<Int, Real>::operator=
    (const MMRGLattice<Int, Real> & lat)
    {
      if (this == &lat)
        return *this;
      this->m_dim = lat.m_dim;
      this->copyBasis(lat);
      this->m_order = lat.m_order;
      //m_ip = lat.m_ip;
      //m_shift = lat.m_shift;
      return *this;
      //MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
      //copy (lat);
      //return *this;
    }

  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::init()
    {
      //kill(); //PW_TODO : wzf ? M-A : Indeed wzf...?
      // This should not be needed
      //m_xi.SetLength(this->m_order);
      //m_A.SetDims(this->m_order, this->m_order);
      //if (this->m_order > ORDERMAX) {
        //m_ip = new bool[1];
        //m_sta.SetDims(1, 1);
      //} else {
        //m_ip = new bool[this->m_order];
        //m_sta.SetDims(this->m_order, this->m_order);
      //}
      //int rmax = std::max(this->m_order, this->getDim());
      //this->m_wSI.SetDims(rmax, this->getDim());

      double temp;
      NTL::conv(temp, this->m_modulo);
      double lgm2 = 2.0 * LatticeTester::Lg (temp);
      //this->calcLgVolDual2 (lgm2);
      //if (m_latType == ORBIT)
      //   initOrbit();
      //PW_TODO
    }

  //===========================================================================

  template<typename Int, typename Real>
    Int & MMRGLattice<Int, Real>::getLac (int j)
    {
      if (isLacunary() && j <= m_lac.getSize() && j > 0)
        return m_lac.getLac(j);
      throw std::out_of_range("MMRGLattice::getLac");
    }


  //===========================================================================
  
  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::setLac(const LatticeTester::Lacunary<Int> & lac)
    {
      m_lac = lac;
      m_lacunaryFlag = true;
      m_numberLacIndices = m_lac.getSize();
    }


  //===========================================================================
  
  //PW_TODO merdasse à tester

  template<typename Int, typename Real>
    std::string MMRGLattice<Int, Real>::toString () const
    {
      std::ostringstream out;
      out << "A = " << m_A << "\n";

      return out.str ();
    }

 //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::buildBasis (int64_t d)
    {
      if (m_lacunaryFlag)
        buildLacunaryBasis(d);
      else
        buildNonLacunaryBasis(d);
    }


  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::buildNonLacunaryBasis (int d)
    // a basis is built in dimension d
    // up to now: implementation only for B = I
    {
      int64_t k, no_powers;
      IntMat A_power;
      k = this->m_A.NumRows();      
      this->m_basis.resize(d, d);
      no_powers = floor(d/k);
      
      
      // Fill the identity matrix in the upper left corner
      for (int64_t i = 0; i < min(k, d); i++) {
         for (int64_t j = 0; j < min(k, d); j++) {
             this->m_basis[i][j] = (i == j);
         }  
      }
      
      // Fill in the powers of A^t which are neede completely 
      A_power = this->m_A; 
      for (int64_t m = 1; m < no_powers; m++) {
          for (int64_t i = 0; i < k; i++) {
               for (int64_t j = 0; j < k; j++) {
                   this->m_basis[i][j+m*k] = A_power[j][i];
               }
          }
          A_power = A_power * this->m_A;
      }
      
      // Fill in the remaining part of (A^{n-1})^t up to dimension d
      for (int64_t i = 0; i < k; i++) {
         for (int64_t j = 0; j < d - no_powers*k; j++) {
             this->m_basis[i][j+no_powers*k] = A_power[j][i];
         }
      }
      
      // Fill m times identity matrix for the rest of the basis matrix.
      for (int64_t i = k; i < d; i++) {
         for (int64_t j = 0; j < d; j++) {
             this->m_basis[i][j] = (i == j) * this->m_modulo;
         }
      }
      
    }

  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::getSubLine(IntVec & vec, IntMat& B, int lign, int jMin,
        int jMax)
    {
      // both jMin and jMax are included
      vec.resize(jMax-jMin+1);
      for (int i = 0; i < (jMax-jMin+1); i++)
        vec[i] = B[lign][jMin+i];
    }

  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::buildLacunaryBasis (int dimension)
    {
    /*
      int sizeA = this->getOrder();
      this->m_vecNorm.resize(dimension);
      this->m_dualvecNorm.resize(dimension);

      int maxIndiceLac = NTL::conv<int>(m_lac[m_lac.getSize()-1]);

      // building the complete basis until: dimension = max lacunary indice
      //-----------------------------------------------------------------------
      IntMat tempBasis;
      //PW_TODO
      tempBasis.resize(this->m_order, maxIndiceLac+1);
      // PW_TODO meilleurs size à trouver
      // +1 because lacunary indices start at 0

      // filling in the diagonal of tempBasis
      for (int i = 0; i < sizeA; i++)
        tempBasis[i][i] = 1;
      //for (int i = sizeA; i < maxIndiceLac; i++)
      //   tempBasis[i][i] = this->m_modulo;

      // using genrator matrix A to complete the first lines of tempBasis
      // with values generated by the recurrence
      NTL::ZZ_p::init(NTL::ZZ(this->m_modulo));
      /
      typename ModInt<Int>::IntMatP temp;
      temp.SetDims(sizeA, sizeA);
      for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

      int maxIter = floor((maxIndiceLac+1)/sizeA);

      for (int k = 1; k < maxIter+1; k++) {
        // calculation of transpose(A^k)
        temp *= NTL::conv<typename ModInt<Int>::IntMatP>(NTL::transpose(m_A));

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
      this->m_wSI.resize(std::max(this->m_order, m_numberLacIndices), m_numberLacIndices);
      this->m_vSI.resize(std::max(this->m_order, m_numberLacIndices), m_numberLacIndices);
      // PW_TODO meilleurs size à trouver

      for (int j = 0; j < m_numberLacIndices; j++) {
        for (int i = 0; i < this->m_order; i++)
          this->m_wSI[i][j] = tempBasis[ i ][ NTL::conv<int>(m_lac[j]) ];
      }
     */
      /*
         std::cout << "tempBasis = \n" << tempBasis << std::endl;
         std::cout << "projection = \n" << this->m_wSI << std::endl;
         std::cout << "\n******************************************\n" << std::endl;
         std::cout << "\nAVANT TRIANGULARIZATION" << std::endl;
         std::cout << "m_vSI = \n" << this->m_vSI << std::endl;
         std::cout << "m_wSI = \n" << this->m_wSI << std::endl;
         */

      // transforming this generating familly into a basis of the lattice
      //-----------------------------------------------------------------------
    /*  LatticeTester::Triangularization <IntMat> (this->m_wSI, this->m_vSI, this->m_order,
          m_numberLacIndices, this->m_modulo);
    */
      //std::cout << "\nAPRES TRIANGULARIZATION, AVANT CALCDUAL" << std::endl;
      //std::cout << "m_vSI = \n" << this->m_vSI << std::endl;
      //std::cout << "m_wSI = \n" << this->m_wSI << std::endl;

      // LatticeTester::CalcDual <IntMat> (this->m_vSI, this->m_wSI, m_numberLacIndices,
       //   this->m_modulo); CW

      /*
         std::cout << "\nAPRES CALCDUAL" << std::endl;
         std::cout << "m_vSI = \n" << this->m_vSI << std::endl;
         std::cout << "m_wSI = \n" << this->m_wSI << std::endl;
         std::cout << "\n******************************************\n" << std::endl;
         */

      //building the basis in dimension 1
    /*  this->m_basis[0][0] = this->m_vSI[0][0];
      this->m_dualbasis[0][0] = this->m_wSI[0][0];
      this->setDim (1);

      this->setNegativeNorm();
      this->setDualNegativeNorm();


      for (int i = 1; i < dimension; i++)
        incDim ();
      */

      // if (!this->checkDuality())
      //   LatticeTester::MyExit (1, "BUG in MMRGLattice::buildNonLacunaryBasis");

    }


  //===========================================================================
  
    template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::incDim()
    {
      if (m_lacunaryFlag)
        incrementDimLacunaryBasis(m_numberLacIndices);
      else
        incrementDimNonLacunaryBasis();
    }

  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::incrementDimNonLacunaryBasis()
    // X_n = A X_{n-1} mod m. We have: dimension >= order.
    {
      IntMat tempdual(this->m_dualbasis);
      LatticeTester::IntLatticeExt<Int, Real>::incDim();
      int newDimension = this->getDim();
      int sizeA = this->getOrder();

      // ************* update of the primal lattice *************
      //  - we add a new coordinate to each vector v_i, this value being
      //  determined by the MMRG recurrence (even if the original vectors have
      //  been transformed linearly and we must apply the same transformations
      //  to their last coordinates).
      //  - we add an extra vector (0,..., 0, m) to complete this dimension
      //    increased basis.


      if (newDimension > this->m_basis_max.NumRows()) {
        // Expanding m_basis_max
        // Taking the last power of m_A computed
        IntMat temp(sizeA, sizeA);
        int start = this->m_basis_max.NumRows();
        for (int i = 0; i < sizeA; i++) {
          for (int j = 0; j < sizeA; j++) {
            temp[i][j] = this->m_basis_max[j][start+i-sizeA];
          }
        }
        // Computing next power of A
        temp = temp * this->m_A;
        IntMat temp2(this->m_basis_max);
        this->m_basis_max.resize(start+sizeA, start+sizeA);
        // Copying old values in new matrix
        for (int i = 0; i < start; i++) {
          for (int j = 0; j < start; j++) {
            this->m_basis_max[i][j] = temp2[i][j];
          }
        }
        // Filling new elements in the matrix
        for (int i = 0; i < sizeA; i++) {
          for (int j = start; j < start+sizeA; j++) {
            this->m_basis_max[i][j] = temp[j-start][i]%this->m_modulo;
          }
          this->m_basis_max[i+start][i+start] = this->m_modulo;
        }
      }

      // Compying the primal
      for (int i = 0; i < newDimension; i++) {
        for (int j = 0; j < newDimension; j++) {
          this->m_basis[i][j] = this->m_basis_max[i][j];
        }
      }
      // Copying the old dual elements
      for (int i = 0; i < newDimension-1; i++) {
        for (int j = 0; j < newDimension-1; j++) {
          this->m_dualbasis[i][j] = tempdual[i][j];
        }
      }

      // Adding the new column and the new row to the dual lattice
      for (int i = 0; i < newDimension-1; i++) {
        this->m_dualbasis[newDimension-1][i] = -this->m_basis[i][newDimension-1];
      }
      this->m_dualbasis[newDimension-1][newDimension-1] = this->m_modulo/this->m_basis[newDimension-1][newDimension-1];


      // if (newDimension <= sizeA) {
      //   this->m_basis[newDimension-1][newDimension-1] = 1;
      // } else {
      //   // we compute the number of steps required to reach the current state of
      //   // the generator for the considered dimension.
      //   int n = floor((newDimension-1) / sizeA);
      //   NTL::ZZ_p::init(NTL::ZZ(this->m_modulo));
      //   typename ModInt<Int>::IntMatP temp;
      //   temp.SetDims(sizeA, sizeA);
      //   for (int i = 0; i < sizeA; i++)
      //     temp[i][i] = 1;
      //   for (int k = 0; k < n; k++)
      //     temp *= NTL::conv<typename ModInt<Int>::IntMatP>(NTL::transpose(m_A));
      //   // PW_TODO : could be useful to keep A^k in memory to shorten computation

      //   // update of the new v_i coordinates using the *temp* matrix and the
      //   // first coefficients of each line (can be seen as a seed vector). So
      //   // this *temp* matrix multiplied by this seed vector gives us the next
      //   // values generated by the MMRG for the considered dimension.
      //   IntVec initialState;
      //   for (int i = 0; i < (newDimension-1); i++) {
      //     getSubLine(initialState, this->m_basis, i, 0, sizeA-1);
      //     initialState = NTL::conv<IntVec>(
      //         NTL::transpose(temp) * NTL::conv<typename ModInt<Int>::IntVecP>(initialState));
      //     this->m_basis[i][newDimension-1] = initialState[newDimension - n*sizeA -1];
      //   }
      //   this->m_basis[newDimension-1][newDimension-1] = this->m_modulo;
      // }

      // // ************* update of the dual basis *************
      // //  - we add a new 0 coordinate to each vector w_i in the dual basis.
      // //  - for the new last line of the matrix, we add an extra vector,
      // //    as described in L'Ecuyer's paper
      // // PW_TODO : citer référence "Guide LatTester"

      // if (newDimension <= sizeA){
      //   this->m_dualbasis[newDimension-1][newDimension-1] = this->m_modulo;
      // } else {

      //   IntVec lastLine;
      //   lastLine.SetLength(newDimension);

      //   for (int i = 0; i < (newDimension-1); i++) {
      //     NTL::matrix_row<const IntMat> row(this->m_dualbasis, i);
      //     lastLine -= this->m_basis[i][newDimension-1] * row;
      //   }
      //   for (int i = 0; i < (newDimension-1); i++)
      //     this->m_dualbasis[newDimension-1][i] = lastLine[i] / this->m_modulo;
      //   this->m_dualbasis[newDimension-1][newDimension-1] = 1;
      // }

      // this->setNegativeNorm();
      // this->setDualNegativeNorm();

      // if (!this->checkDuality())
      //   LatticeTester::MyExit (1, "BUG in MMRGLattice::incrementDimBasis");

    }

  //===========================================================================

  template<typename Int, typename Real>
    void MMRGLattice<Int, Real>::incrementDimLacunaryBasis(int Imax)
    {
      LatticeTester::IntLatticeExt<Int, Real>::incDim();
      const int dim = this->getDim (); // new dimension (dim++)

      //PW_TODO
      /*
         if (dim >= IMax) {
         MyExit (0,
         "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
         }
         */

      IntVec tempLineBasis (dim);
      IntVec tempColBasis (dim);

      for (int i = 0; i < dim-1; i++) {

        // tempLineBasis <- m_basis[i]
        for (int k = 0; k < dim-1; k++)
          tempLineBasis[k] = this->m_basis[i][k];

        // v[i] -> VSI[0].
        // for (int j = 0; j < dim; j++)
        //    this->m_vSI[0][j] = this->m_basis[i][j];

        //clear (this->m_vSI[i][0]);

        for (int i1 = 0; i1 < dim-1; i1++) {

          Int tempScalDual;

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

        //clear (this->m_wSI[0][j]);
        Int tempScalDualBis;

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

      //this->setDim (dim + 1);
      this->setNegativeNorm ();
      this->setDualNegativeNorm ();

    }


    
 
} // End namespace LatMRG
#endif

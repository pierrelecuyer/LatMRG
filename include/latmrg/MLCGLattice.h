#ifndef LATMRG_MLCGLATTICE_H
#define LATMRG_MLCGLATTICE_H

#include <cassert>
//#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FlexTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {

/**
 * This class is similar to `LCGLattice` and `MRGLattice`, but for matrix LCGs (MLCGs).
 * An MLCG is defined by selecting a \f$k\times k\f$ generator matrix \f$A\f$ and
 * a \f$w\times k\f$ output matrix \f$B\f$.
 * Its state is a \f$k\f$-bit vector \f$\mathbb{x}\f$ that follows the recurrence
 * \f[
 *   \mathbb{x}_n = \mathbb{A} \mathbb{x}_{n-1} \mod m
 * \f]
 * and the output is a \f$w\f$-bit vector \f$\mathbb{y}\f$ defined by
 * \f[
 *   \mathbb{y}_n = \mathbb{B} \mathbb{x}_{n} \mod m.
 * \f]
 * When \f$B\f$ is not given, it is assumed to be the identity and this
 * simplifies most of the functions.
 * There are constructors that do not take \f$A\f$ and/or \f$B\f$ as inputs,
 * but only their sizes, so the same `MLCGLattice` object can be used for
 * several different choices of matrices, set by `setAB` or `setA`.
 * Note that an MRG is a special case of a MLCG for which \f$A\f$ has a
 * special form, \f$w=1\f$, and \f$B\f$ is just a unit row vector.
 *
 * *The current implementation works only for \f$B = I\f$.*
 */
template<typename Int, typename Real>
class MLCGLattice: public LatticeTester::IntLatticeExt<Int, Real> {

public:

   /**
    * This constructor takes as input the modulus `m`, the matrices `A` and `B`,
    * some maximal dimensions, and the norm used to measure the vector lengths.
    * The values of `k` and `w` are deduced from the numbers of rows of `A` and `B` and cannot be changed afterwards.
    * The parameters `maxDim` and  `maxDimProj` give the maximal dimension of a basis (primal or m-dual) and
    * the maximal dimension of a projection (for `buildProjection`).
    * The `maxDim` should be a multiple of `w`.
    * The constructor does not build a basis.
    */
   MLCGLattice(const Int &m, const IntMat &A, const IntMat &B, int64_t maxDim, int64_t maxDimProj,
         NormType norm = L2NORM);

   /**
    * This version assumes that `B` is the `k x k` identity.
    */
   MLCGLattice(const Int &m, const IntMat &A, int64_t maxDim, int64_t maxDimProj, NormType norm =
         L2NORM);

   /**
    * Same as the first constructor, except that only the dimensions `k` and `w` are given instead of the
    * matrices `A` and `B`.  These matrices can be set via `setAB` of `setA`.
    */
   MLCGLattice(const Int &m, int64_t k, int64_t w, int64_t maxDim, int64_t maxDimProj,
         NormType norm = L2NORM);

   /**
    * Same as the previous constructor, except that `B` is assumed to be the identity.
    */
   MLCGLattice(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj, NormType norm = L2NORM);

   /**
    * Destructor.
    */
   virtual ~MLCGLattice();

   /**
    * Cleans and releases memory used by this object.
    */
   virtual void kill();

   /**
    * Sets the matrices `A` and `B`. They must be of the right sizes.
    */
   // virtual void setAB(const IntMat &A, const IntMat &B);

   /**
    * Sets only the matrix `A` and assumes that `B` is the `k x k` identity matrix.
    */
   virtual void setA(const IntMat &A);

   /**
    * Returns a non-mutable copy of the generator matrix `A`.
    */
   const IntMat& getA() const {
      return m_A;
   }

   /**
    * Builds a primal basis in `dim` dimensions, where `dim` must be at least `k` and
    * not exceed `m_maxDim`.
    * This initial primal basis will be upper triangular.
    * This function actually builds an original basis in up to `maxDim` dimensions,
    * and initializes the current basis to `dim` dimensions.
    */
   virtual void buildBasis(int64_t dim);

   /**
    * Builds a lower-triangular m-dual basis directly in `dim` dimensions.
    * This `dim` must be at least `k`and  must not exceed `m_maxDim`.
    * This function actually builds an original dual basis in up to `maxDim` dimensions,
    * and initializes the current dual basis to `dim` dimensions.
    */
   virtual void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `m_maxDim`.
    */
   virtual void incDimBasis();

   /**
    * Increases the current dimension of the m-dual basis by 1.
    * The new increased dimension must not exceed `m_maxDim`.
    */
   virtual void incDimDualBasis();

   /**
    * Builds a basis for the projection of this `MLCGLattice` onto the coordinates
    * in `coordSet` and puts it as the `m_basis` of `projLattice`.
    * The largest coordinate of the projection must not exceed `m_maxDim` and
    * the number of coordinates must not exceed `m_maxDimProj`.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` in `projLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    *  Compute and store the powers of `A` mod m, from 0 to `numPow-1`,
    *  and put the transposed powers in a matrix.
    *  These powers are used to build a basis or m-dual basis.
    *  When some powers have been computed already, the function computes only the missing ones.
    */
   void computePowersOfA(int64_t numPow);

protected:

   /**
    * True iff `B` is not the identity.
    */
   bool m_withB = false;

   /**
    * The number of rows in `A` and in `B`.
    */
   int64_t m_k, m_w;

   /**
    * The generator and output matrices of the recurrence.
    */
   IntMat m_A, m_B;

   /**
    * These objects store A^p and B A^p,
    * where p = m_pow is the highest power of A computed so far.
    * This power is updated when we call `computePowersOfA`.
    */
   typename FlexModInt<Int>::IntMatP m_Am;    // Matrix A in Z_m.
   typename FlexModInt<Int>::IntMatP m_powA;  // Matrix A^p in Z_m.
   int64_t m_pow = 0;
   // IntMat m_BpowA;
   // IntMat m_tempA;

   /**
    * Maximal number of coordinates allowed in a projection.
    */
   int64_t m_maxDimProj;

   /**
    * Used to store a copy of the original matrix basis and/or its m-dual.
    * These matrices are initialized by `buildBasis` or `duildDualBasis`,
    * and expanded when needed.
    */
   IntMat m_powersOfA;
   // IntMat m_dualbasisOriginal;

   /**
    * This auxiliary matrix is used to store the set of generating vectors of a projection
    * before reducing them into a triangular basis.
    * Its size must be at least `(d + m_k) x d`, where `d = maxDimProj` is the largest dimension of a projection.
    * It is used only in `buildProjection`.
    */
   IntMat m_genTemp;

   /**
    * Indicates which lattice or sublattice is analyzed.
    */
   // LatticeType m_latType;
};
// End class declaration

//===========================================================================
// IMPLEMENTTION

//=============================================================================

// Main constructor.
template<typename Int, typename Real>
MLCGLattice<Int, Real>::MLCGLattice(const Int &m, int64_t k, int64_t w, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   m_k = k;
   if (w == 0) {
      m_withB = false;
      m_w = k;  // B is the identity.
   } else {
      m_withB = true;
      m_w = w;
   }
   FlexModInt<Int>::mod_init(this->m_modulo);  // Set `m` for modular arithmetic.
   assert(((maxDim % m_w) == 0) && "maxDim must be a multiple of m_w");
   m_maxDimProj = maxDimProj;
   m_genTemp.SetDims(maxDimProj + k, maxDimProj);
   m_pow = 0;
   // m_powA.SetDims(k, k);
   // m_tempA.SetDims(k, k);
   m_powersOfA.SetDims(k, maxDim);
}

template<typename Int, typename Real>
MLCGLattice<Int, Real>::MLCGLattice(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj,
      NormType norm) :
      MLCGLattice<Int, Real>(m, k, 0, maxDim, maxDimProj, norm) {
}

template<typename Int, typename Real>
MLCGLattice<Int, Real>::MLCGLattice(const Int &m, const IntMat &A, const IntMat &B, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLattice<Int, Real>(m, A.NumRows(), B.NumRows(), maxDim, maxDimProj) {
   // setAB(A, B);
   assert((false) && "The constructor with a matrix `B` is not implemented");
}

template<typename Int, typename Real>
MLCGLattice<Int, Real>::MLCGLattice(const Int &m, const IntMat &A, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLattice<Int, Real>(m, A.NumRows(), 0, maxDim, maxDimProj, norm) {
   setA(A);
}

//===========================================================================

template<typename Int, typename Real>
MLCGLattice<Int, Real>::~MLCGLattice() {
   kill();
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::kill() {
   m_A.kill();
   m_B.kill();
   m_Am.kill();
   m_powA.kill();
   m_powersOfA.kill();
   m_genTemp.kill();
   // m_bV0.kill();
}

//============================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::setA(const IntMat &A) {
   assert(A.NumRows() == this->m_k);  // `A` must have the right size.
   this->m_A = A;
   m_Am = conv<typename FlexModInt<Int>::IntMatP>(m_A);
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   m_pow = 0;
}

//============================================================================
/*
template<typename Int, typename Real>
void MLCGLattice<Int, Real>::setAB(const IntMat &A, const IntMat &B) {
   setA(A);
   assert(B.NumCols() == this->m_k);  // `B` must have the right size.
   assert(B.NumRows() == this->m_w);
   this->m_B = B;
}
*/
//===========================================================================
// Computes the first `k` rows of the initial basis in up to `numPow * k` dimensions,
// under the assumption that B = I.  This is (I, Y_1^\tr, ..., Y_{n-1}^\tr).
// This matrix is used by `buildBasis`, `buildBasisDual`, `incDimDualBasis`,
// but not by`incDimBasis`.
template<typename Int, typename Real>
void MLCGLattice<Int, Real>::computePowersOfA(int64_t numPow) {
   assert(!m_withB);
   assert(numPow > 0);                      // At least w dimensions.
   assert(numPow * m_w <= this->m_maxDim);  // No more than maxDim dimensions.
   int64_t i, j;
   if (numPow <= m_pow) return;  // Nothing to do.
   if (m_pow == 0) {
      NTL::ident(m_powA, m_k);   // For first block, must be identity matrix = A^0.
      // if (m_withB) pBA = m_B;     // There is a B matrix.
      // Put identity in the upper left block.
      for (i = 0; i < m_k; i++) {
         for (j = 0; j < m_k; j++) {
            m_powersOfA[i][j] = (i == j);
         }
      }
      m_pow = 1;
   }
   // Compute the blocks in the first row, in succession.
   for (int64_t p = m_pow; p < numPow; p++) {
/*
       for (i = 0; i < m_k; i++) {
         for (j = 0; j < m_k; j++) {
            m_tempA[i][j] = 0;
            for (int64_t l = 0; l < m_k; l++)
               m_tempA[i][j] += m_A[i][l] * m_powA[l][j];
            m_tempA[i][j] = m_tempA[i][j] % this->m_modulo;
         }
      }
      m_powA = m_tempA;
*/
      m_powA = m_powA * m_Am;
      for (i = 0; i < m_k; i++)
         for (j = 0; j < m_w; j++)
            m_powersOfA[i][p * m_k + j] = conv<Int>(m_powA[j][i]); // transposed.
   }
   m_pow = numPow;
}
/*
//===========================================================================
// Computes the first `k` rows of the initial basis in up to `numPow * k` dimensions,
// under the assumption that B = I.  This is (I, Y_1^\tr, ..., Y_{n-1}^\tr).
// This matrix is used by `buildBasis`, `buildBasisDual`, `incDimDualBasis`,
// but not by`incDimBasis`.
template<typename Int, typename Real>
void MLCGLattice<Int, Real>::computePowersOfA(int64_t numPow) {
   assert(!m_withB);
   assert(numPow > 0);                      // At least w dimensions.
   assert(numPow * m_w <= this->m_maxDim);  // No more than maxDim dimensions.
   int64_t i, j;
   if (numPow <= m_pow) return;  // Nothing to do.
   if (m_pow == 0) {
      NTL::ident(m_powA, m_k);   // For first block, must be identity matrix = A^0.
      // if (m_withB) pBA = m_B;     // There is a B matrix.
      // Put identity in the upper left block.
      for (i = 0; i < m_k; i++) {
         for (j = 0; j < m_k; j++) {
            m_powersOfA[i][j] = (i == j);
         }
      }
      m_pow = 1;
   }
   // Compute the blocks in the first row, in succession.
   for (int64_t p = m_pow; p < numPow; p++) {
      for (i = 0; i < m_k; i++) {
         for (j = 0; j < m_k; j++) {
            m_tempA[i][j] = 0;
            for (int64_t l = 0; l < m_k; l++)
               m_tempA[i][j] += m_A[i][l] * m_powA[l][j];
            m_tempA[i][j] = m_tempA[i][j] % this->m_modulo;
         }
      }
      m_powA = m_tempA;
      for (i = 0; i < m_k; i++)
         for (j = 0; j < m_w; j++)
            m_powersOfA[i][p * m_k + j] = m_powA[j][i]; // transposed.
   }
   m_pow = numPow;
}
*/

//===========================================================================

// Up to now: implementation works only for B = I.
template<typename Int, typename Real>
void MLCGLattice<Int, Real>::buildBasis(int64_t dim) {
   assert(!m_withB);
   this->setDim(dim);
   int64_t numPow = dim / m_w;
   if (m_w * numPow < dim) numPow++;
   int i, j, l;
   computePowersOfA(numPow);
   for (int64_t i = 0; i < m_k; i++)
      for (int64_t j = 0; j < dim; j++) {
         if (!m_withB) this->m_basis[i][j] = m_powersOfA[i][j];  // No B.
         else {  // For the future...
            // m_BpowA = m_B * m_powA % this->m_modulo;
            Int sum = Int(0);
            for (l = 0; l < m_k; l++)
               sum += m_B[i][l] * m_powersOfA[i][l + (j / m_w) * m_k];
            this->m_basis[i][j] = sum % this->m_modulo;
         }
      }
   // Fill m times the identity matrix for the rest of the basis.
   // This assumes that B=I.
   for (i = m_k; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         this->m_basis[i][j] = (i == j) * this->m_modulo;
      }
   }
   // std::cout << "buildBasis, basis = \n" << this->m_basis << "\n";
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::incDimBasis() {
   assert(!m_withB);
   int64_t d = this->getDim();      // Current dimension.
   int64_t p = d / m_k - 1;         // Previous block number, starts at 0.
   int64_t r = d - m_k * (p + 1);   // Position in block, starts  at 0.
   Int sum;
   // Add a new coordinate to each basis vector.
   for (int64_t i = 0; i < d; i++) {
      sum = Int(0);
      for (int64_t j = 0; j < m_k; j++)
         sum += this->m_basis[i][p * m_k + j] * m_A[r][j];
      this->m_basis[i][d] = sum;
   }
   // Then add a new row vector.
   for (int64_t j = 0; j < d; j++) {
      this->m_basis[d][j] = 0;
   }
   this->m_basis[d][d] = this->m_modulo;
   // std::cout << "m_basis after incDimBasis = \n" << this->m_basis << "\n";
   this->setDim(d + 1);
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::buildDualBasis(int64_t dim) {
   assert(!m_withB);
   computePowersOfA(this->m_maxDim / m_w);
   this->setDimDual(dim);
   for (int64_t i = 0; i < m_k; i++)
      for (int64_t j = 0; j < dim; j++)
         this->m_dualbasis[i][j] = (i == j) * this->m_modulo;
   for (int64_t i = m_k; i < dim; i++) {
      for (int64_t j = 0; j < m_k; j++)
         this->m_dualbasis[i][j] = -m_powersOfA[j][i];
      for (int64_t j = m_k; j < dim; j++)
         this->m_dualbasis[i][j] = (i == j);
   }
   // std::cout << "buildDualBasis, dualbasis = \n" << this->m_dualbasis << "\n";
}

//===========================================================================

// This function assumes that `buildBasisOriginal` has been computed before,
// in up to at least the new dimension.
template<typename Int, typename Real>
void MLCGLattice<Int, Real>::incDimDualBasis() {
   assert(!m_withB);
   int64_t d = this->getDimDual();  // New current dimension.
   // Add a new zero coordinate to each dual basis vector.
   for (int64_t i = 0; i < d; i++)
      this->m_dualbasis[i][d] = 0;
   this->m_dualbasis[d][d] = 1;
   for (int64_t j = m_k; j < d; j++)
      this->m_dualbasis[d][j] = 0;
   for (int64_t j = 0; j < m_k; j++)
      this->m_dualbasis[d][j] = -m_powersOfA[j][d];
   // std::cout << "m_basis after incDimDualBasis = \n" << this->m_dualbasis << "\n";
   this->setDimDual(d + 1);
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   assert(!m_withB);   // Implemented only for B=I for now.
   int64_t d = coordSet.size();
   assert(d <= m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(this->m_maxDim));
   projLattice.setDim(d);
   IntMat* pbasis = &projLattice.getBasis(); // Pointer to basis of projection.
   int64_t i, j;
   int64_t k = this->m_k;
   bool projCase1 = true; // This holds if the first m_k coordinates are all in `coordSet`.
   // Check if we are in case 1.
   // This assumes that the coordinates of each projection are always in increasing order!  ***
   if (d < (unsigned) k) projCase1 = false;
   else {
      j = 1;
      for (auto it = coordSet.begin(); it != coordSet.end(), j <= k; it++, j++) {
         if (*it != unsigned(j)) {
            projCase1 = false;
            break;
         }
      }
   }
   int64_t iadd = 0;
   if (!projCase1) {
      iadd = k;  // We will have k more rows in pbasis.
      pbasis = &m_genTemp;  // Pointer will point a set of gen vectors.
   }
   // We first copy the selected coordinates for first k rows.
   j = 0;
   for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++)
      for (i = 0; i < k; i++)
         (*pbasis)[i][j] = m_powersOfA[i][*it - 1];
   // Then the other rows.
   for (i = k; i < d + iadd; i++)
      for (j = 0; j < d; j++)
         (*pbasis)[i][j] = this->m_modulo * (i == j + iadd);
   // If not case1, we must reduce the set of gen vectors.
   if (!projCase1)
      upperTriangularBasis(projLattice.getBasis(), m_genTemp, this->m_modulo, this->m_dim, d);
}

//============================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   assert(!m_withB);
   int64_t d = coordSet.size();     // The dimension of this projection.
   assert(d <= m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(this->m_maxDim));
   projLattice.setDimDual(d);
   // We build a basis for the projection and we compute its m-dual.
   this->buildProjection(projLattice, coordSet);
   mDualUpperTriangular(projLattice.getDualBasis(), projLattice.getBasis(), this->m_modulo, d);
}

//===========================================================================

}// End namespace LatMRG
#endif

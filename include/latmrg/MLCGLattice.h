#ifndef LATMRG_MLCGLATTICE_H
#define LATMRG_MLCGLATTICE_H

#include <cassert>
//#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FlexTypes.h"
#include "latticetester/IntLatticeExt.h"
//#include "latticetester/MRGLattice.h"
//#include "latmrg/FlexModInt.h"
//#include "Primitivity.h"

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
 * several different choices of matrices, set by `setAB`.
 * Note that an MRG is a special case of a MLCG for which \f$A\f$ has a
 * special form, \f$w=1\f$, and \f$B\f$ is just a unit row vector.
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
   ~MLCGLattice();

   /**
    * Cleans and releases memory used by this object.
    */
   void kill();

   /**
    * Sets the matrices `A` and `B`. They must be of the right sizes.
    */
   void setAB(const IntMat &A, const IntMat &B);

   /**
    * Sets only the matrix `A` and assumes that `B` is the `k x k` identity matrix.
    */
   void setA(const IntMat &A);

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
    * in `proj` and puts it as the `m_basis` of `projLattice`.
    * The largest coordinate of the projection must not exceed `m_maxDim` and
    * the number of coordinates must not exceed `m_maxDimProj`.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` in `projLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;


protected:

   /**
    *  Compute and store the powers of `A` mod m, from 0 to `numPow-1`,
    *  and put the transposed powers in a matrix.
    *  These powers are used to build a basis or m-dual basis.
    */
   void computePowersOfA(int64_t numPow);

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
    * These objects store A^p, B A^p, and a temporary matrix to compute them,
    * where p = m_pow is the highest power of A computed so far
    */
   IntMat m_powA;
   // IntMat m_BpowA;
   IntMat m_tempA;
   int64_t m_pow = 0;

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
   std::cout << "maxDim = " << maxDim << ", k = " << k << ", m_w = " << m_w << "\n";
   assert(((maxDim % m_w) == 0) && "maxDim must be a multiple of m_w");
   m_maxDimProj = maxDimProj;
   m_genTemp.SetDims(maxDimProj + k, maxDimProj);
   m_pow = 0;
   m_powA.SetDims(k, k);
   m_tempA.SetDims(k, k);
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
   m_A = A;
   m_B = B;
}

template<typename Int, typename Real>
MLCGLattice<Int, Real>::MLCGLattice(const Int &m, const IntMat &A, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLattice<Int, Real>(m, A.NumRows(), 0, maxDim, maxDimProj, norm) {
   m_A = A;
}

//===========================================================================

template<typename Int, typename Real>
MLCGLattice<Int, Real>::~MLCGLattice() {
   kill();
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::kill() {
// To do...
}

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
   m_pow = 0;
   NTL::ident(m_powA, m_k);   // For first block, must be identity matrix = A^0.
   // IntMat & pA = m_powA;      // Reference to a matrix, either m_powA or m_BpowA.
   // if (m_withB) pBA = m_B;     // There is a B matrix.

   // Put identity in the upper left block.
   for (i = 0; i < m_k; i++) {
      for (j = 0; j < m_k; j++) {
         m_powersOfA[i][j] = (i == j);
      }
   }
   // Compute the blocks in the first row, in succession.
   for (int64_t p = 1; p < numPow; p++) {
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

//===========================================================================

// up to now: implementation works only for B = I.
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
               sum += m_B[i][l] * m_powersOfA[i][l + (j/m_w)*m_k];
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
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::incDimBasis() {
   assert(!m_withB);
   int64_t d = this->getDim();      // Current dimension.
   int64_t p = d / m_k - 1;         // Previous block number, starts at 0.
   int64_t r = d - m_k * (p + 1);   // Position in block, starts  at 0.
   Int sum;
   // Add a new zero coordinate to each basis vector.
   for (int64_t i = 0; i < d; i++) {
      sum = Int(0);
      for (int64_t j = 0; j < m_k; j++)
         sum += this->m_basis[i][p * m_k + j] * m_A[r][j];
      this->m_basis[i][d] = sum;
   }
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
      const Coordinates &proj) {
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
}

}   // End namespace LatMRG
#endif

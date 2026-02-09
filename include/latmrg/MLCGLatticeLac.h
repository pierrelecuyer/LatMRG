#ifndef LATMRG_MLCGLATTICELAC_H
#define LATMRG_MLCGLATTICELAC_H

#include "latticetester/EnumTypes.h"
// #include "latticetester/Types.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/MLCGLattice.h"
#include "latmrg/FlexModInt.h"
#include "Primitivity.h"

namespace LatMRG {

/**
 * This subclass of `MLCGLattice` constructs and handles lattice with arbitrary lacunary
 * indices that can be spaced far apart, similar to `MRGLatticeLac`.
 * To construct a basis, we compute the \f$k\times s\f$ matrix \f$\mathbf{V}_{I,k\times s}\f$
 * as explained in Section 3.2.5 of the guide.
 * For this, we compute the large powers of the \f$k\times k\f$ matrix \f$A\f$ modulo \f$m\f$ directly.
 * When the dimension \f$k\f$ of \f$A\f$ is small, this is good enough.
 * For larger \f$k\f$, it is more efficient to use the polynomial representation
 * of the state, as explained in the guide.
 *
 * *This implementation assumes that the matrix \f$B\f$ is the identity.*
 *
 * ....   TO DO: Implement both and compare.   ....
 */

template<typename Int, typename Real>
class MLCGLatticeLac: public MLCGLattice<Int, Real> {

public:

   /**
    * These constructors have the same parameters as in `MLCGLattice` and operate in a similar way.
    */
   MLCGLatticeLac(const Int &m, const IntMat &A, const IntMat &B, int64_t maxDim,
         int64_t maxDimProj, NormType norm = L2NORM);

   /**
    * This version assumes that `B` is the `k x k` identity.
    */
   MLCGLatticeLac(const Int &m, const IntMat &A, int64_t maxDim, int64_t maxDimProj, NormType norm =
         L2NORM);

   /**
    * Here, only the dimensions `k` and `w` are given instead of the matrices `A` and `B`.
    */
   MLCGLatticeLac(const Int &m, int64_t k, int64_t w, int64_t maxDim, int64_t maxDimProj,
         NormType norm = L2NORM);

   /**
    * Here, `B` is assumed to be the identity.
    */
   MLCGLatticeLac(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj, NormType norm =
         L2NORM);

   /**
    * Destructor.
    */
   virtual ~MLCGLatticeLac();

   /**
    * Cleans and releases memory used by this object.
    */
   virtual void kill();

   /**
    * Sets the matrix `A` and assumes that `B` is the `k x k` identity matrix.
    */
   virtual void setA(const IntMat &A);

   /**
    * Sets the vector of lacunary indices to `lac`, whose length should be at least `k` and
    * not exceed `m_maxDim`. The indices are assumed to be in increasing order.
    */
   void setLac(const IntVec &lac);

   /**
    * Returns the \f$j\f$-th lacunary index,\f$i_j\f$.
    */
   Int& getLac(int j) {
      return m_lac[j];
   }
   ;

   /**
    * Computes the first `dim` columns of the matrix \f$\mathbb{V}_{I,k\times s}\f$
    * defined in the guide, for the current matrices \f$A\f$ and \f$B\f$ and the current
    * set of lacunary indices \f$I = \{i_1,\dots,i_s\}\f$.
    * If some columns were already computed, it computes only the missing ones.
    * When the parameter `dim` is not given (with the default value of 0),
    * the number of columns is set to the number `m_s` of lacunary indices.
    * This function is called internally when needed, when we build a basis.
    * The first version computes powers of `A` directly, the second uses polynomial arithmetic.
    */
   virtual void computeVIks(int64_t dim = 0);

   // void computeVIksPoly(int64_t dim = 0);

   /**
    * Builds an initial upper-triangular primal basis in `dim` dimensions.
    * This `dim` must be at least `k` and must not exceed `m_maxDim`.
    * The function `computeVIks` is called internally if needed.
    */
   virtual void buildBasis(int64_t dim) override;

   /**
    * Builds a lower-triangular m-dual basis in `dim` dimensions.
    * This function first calls `buildBasis` and inverts the basis to get its m-dual.
    */
   virtual void buildDualBasis(int64_t dim) override;

   /**
    * Increases the current dimension of the primal basis by 1 and updates it.
    * The new increased dimension must not exceed `m_maxDim`.
    */
   virtual void incDimBasis() override;

   /**
    * Increases the current dimension of the m-dual basis by 1.
    * The new increased dimension must not exceed `m_maxDim`.
    * This function adds a zero coordinate to each vector of the current m-dual basis,
    * and it adds the last row of the matrix \f$\mathbf{W}_{I,s}\f$ which is the m-dual
    * of the matrix \f$\mathbf{V}_{I}\f$ in \f$s\f$ dimensions.
    */
   virtual void incDimDualBasis() override;

   /**
    * Builds a basis for the projection of the current `MLCGLatticeLac` onto the coordinates
    * in `coordSet` and puts it as the `m_basis` of `projLattice`.
    * This is done simply by projecting the current `m_basis`, as in the default
    * implementation in `IntLattice`.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` in `projLattice`.  Also done as the default in `IntLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

protected:

   /**
    * Initializes some of the local variables.
    */
   // void init();
   /**
    * Builds an initial upper-triangular basis in `dim` dimensions and puts it in `basis`.
    */
   void buildBasisOriginal(IntMat &basis, int64_t dim);

   /**
    * Stores the set \f$I\f$ of lacunary indices and its length \f$s\f$.
    */
   IntVec m_lac;

   int64_t m_s = 0;

   /**
    * This variable will be true iff if the set I contains {1,...,k}.
    * It is set when we set the lacunary indices and used when we
    * compute `m_VIks` or build a basis.
    */
   bool m_case1 = true;

   /**
    * This is the matrix \f$\mathbb{V}_{I,k\times s}\f$ defined in the guide.
    * It contains the ingredients we need to build a basis.
    * Each time the vector of coefficients \f$\mathbb{a}\f$ or the set \f$I\f$ of
    * lacunary indices changes, we must recompute this matrix.
    * This is done by `computeVIks`.
    */
   IntMat m_VIks;

   /**
    * The current number of columns for which `m_VIks` is computed.
    */
   int64_t m_dimVIks = 0;

   /**
    * Auxiliary matrix used to store a set of generating vectors used to compute a
    * triangular basis when `m_case1 == false`.
    * Its size must be at least `(maxDim + k) x maxDim`.
    * It is already defined in the parent class, but its size here may be larger.
    */
   // IntMat m_genTemp;
   /**
    * Used to store a copy of the original matrix basis and its m-dual,
    * in case we examine projections over successive coordinates by adding the
    * coordinates incrementally via `incDimBasis` or `incDimDualBasis`.
    * These matrices are initialized by `buildBasis` and `duildDualBasis`.
    */
   IntMat m_basisOriginal;
   IntMat m_dualbasisOriginal;

   IntMat m_M;  // Matrix M used in `incDimBsis`.
   int64_t m_dimM = 0;  // Its current dimension.

   /**
    * The characteristic polynomial P(z) of `A`, set in `setA`.
    */
   typename FlexModInt<Int>::PolX m_Pz;

   /**
    * Matrices used to compute `m_VIks`.
    */
   //typename FlexModInt<Int>::IntMatP m_Am;   // Matrix A in mod m.
   //typename FlexModInt<Int>::IntMatP m_powA; // Matrix A to some power, in mod m.
   typename FlexModInt<Int>::IntMatP m_Ym;   // Matrix Y in mod m.

};
// End class declaration

//===========================================================================
// Implementation:

//===========================================================================

template<typename Int, typename Real>
MLCGLatticeLac<Int, Real>::MLCGLatticeLac(const Int &m, int64_t k, int64_t w, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLattice<Int, Real>(m, k, w, maxDim, maxDimProj, norm) {
   m_VIks.SetDims(k, maxDim);
   this->m_genTemp.SetDims(maxDim + k, maxDim);
   m_basisOriginal.SetDims(maxDim, maxDim);
   m_dualbasisOriginal.SetDims(maxDim, maxDim);
   m_M.SetDims(maxDim, maxDim);
   FlexModInt<Int>::mod_init(this->m_modulo);  // Set `m` for modular arithmetic.
}

template<typename Int, typename Real>
MLCGLatticeLac<Int, Real>::MLCGLatticeLac(const Int &m, int64_t k, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLatticeLac<Int, Real>(m, k, 0, maxDim, maxDimProj, norm) {
}

template<typename Int, typename Real>
MLCGLatticeLac<Int, Real>::MLCGLatticeLac(const Int &m, const IntMat &A, const IntMat &B,
      int64_t maxDim, int64_t maxDimProj, NormType norm) :
      MLCGLatticeLac<Int, Real>(m, A.NumRows(), B.NumRows(), maxDim, maxDimProj, norm) {
   this->setAB(A, B);
}

template<typename Int, typename Real>
MLCGLatticeLac<Int, Real>::MLCGLatticeLac(const Int &m, const IntMat &A, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MLCGLatticeLac<Int, Real>(m, A.NumRows(), maxDim, maxDimProj, norm) {
   this->setA(A);
}

//===========================================================================

template<typename Int, typename Real>
MLCGLatticeLac<Int, Real>::~MLCGLatticeLac() {
   kill();
}

//===========================================================================

template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::kill() {
}

//============================================================================
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::setLac(const IntVec &lac) {
   int64_t k = this->m_k;
   // std::cout << "\n setLac, k = " << k << ",  maxDim = " << this->m_maxDim << ",\n lac = " << lac << "\n";
   assert((lac.length() >= k) && (lac.length() <= this->m_maxDim));
   m_lac = lac;
   m_s = lac.length();
   m_dimVIks = 0;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   // Check if set I contains {1,...,k} (case1 is true).
   if (lac[k - 1] == k) this->m_case1 = true;
   else this->m_case1 = false;
   // std::cout << "Are we in case 1? " << m_case1 << "\n";
}

//============================================================================

template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::setA(const IntMat &A) {
   assert(A.NumRows() == this->m_k);  // `aa` must have the right length.
   this->m_A = A;
   m_dimVIks = 0;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   // Set the characteristic polynomial P(z) of the recurrence in NTL.
   this->m_Am = conv<typename FlexModInt<Int>::IntMatP>(A);
   getCharacPoly<Int>(m_Pz, this->m_Am);
   FlexModInt<Int>::PolE::init(m_Pz);
}

//============================================================================

// This function computes the first `dim` columns of the matrix `m_VIks`
// by computing the relevant powers of `A`.
// If some columns were already computed, it computes only the missing ones.
// By default, `dim` is set to the number of lacunary indices.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::computeVIks(int64_t dim) {
   assert(dim <= min(this->m_maxDim, m_s));  // Upper bound on `dim`.
   if (dim == 0) dim = m_s;
   if (dim <= m_dimVIks) return;     // Nothing to do.
   int64_t k = this->m_k;
   int64_t i, j, r;
   Int n;
   Int nprev = Int(0);  // Current and previous power n.
   ident(m_Ym, k);
   // std::cout << "First m_Yp = \n" << m_Yp << "\n";
   for (j = m_dimVIks; j < dim; j++) {
      // std::cout << "Inside computeVIks, j = " << j << "\n";
      n = (m_lac[j] - 1) / k;
      r = (m_lac[j] - 1) % k;
      // std::cout << "  n = " << n << ", r = " << r << "\n";
      if (n > nprev) {
         // std::cout << "  n - nprev = " << n - nprev << "\n";
         power(this->m_powA, this->m_Am, n - nprev);
         mul(m_Ym, this->m_powA, m_Ym);
         // std::cout << "New m_Yp = \n" << m_Yp << "\n";
         nprev = n;
      }
      // std::cout << "Inside computeVIks, m_Ym = \n" << m_Ym << "\n";
      for (i = 0; i < k; i++)
         m_VIks[i][j] = conv<Int>(m_Ym[r][i]);
   }
   m_dimVIks = dim;
   // std::cout << "Finished computeVIks: m_VIks = \n" << m_VIks << "\n";
}

//============================================================================

// Builds an upper-triangular `basis` in `dim` dimensions, using `m_VIks`.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::buildBasisOriginal(IntMat &basis, int64_t dim) {
   int64_t k = this->m_k;
   assert(dim >= k);
   assert(dim <= min(this->m_maxDim, m_s));
   if (dim > m_dimVIks) this->computeVIks(dim);
   int64_t i, j;
   // Fill the first k rows using m_VIks.
   IntMat* pbasis = &basis;      // reference to a basis.
   // std::cout << "buildBasisOriginal1, pbasis = \n" << pbasis << "\n";
   if (!m_case1) pbasis = &this->m_genTemp;
   for (i = 0; i < k; i++)
      for (j = 0; j < dim; j++)
         (*pbasis)[i][j] = m_VIks[i][j];
   // Fill the other rows.
   if (m_case1) {
      for (i = k; i < dim; i++)
         for (j = 0; j < dim; j++)
            (*pbasis)[i][j] = (i == j) * this->m_modulo;
   } else {
      for (i = 0; i < dim; i++)
         for (j = 0; j < dim; j++)
            (*pbasis)[i + k][j] = (i == j) * this->m_modulo;
      // std::cout << "buildBasisOriginal2, gen vectors = \n" << pbasis << "\n";
      upperTriangularBasis(m_basisOriginal, *pbasis, this->m_modulo, dim + k, dim);
   }
   //std::cout << "buildBasisOriginal4, basis = \n" << m_basisOriginal << "\n";
}

//============================================================================

// Builds an upper-triangular `basis` in `dim` dimensions, using `m_VIks`.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::buildBasis(int64_t dim) {
   this->setDim(dim);
   // if (m_dimVIks < this->m_maxDim) computeVIks(this->m_maxDim);
   buildBasisOriginal(m_basisOriginal, this->m_maxDim);
   for (int64_t i = 0; i < dim; i++)
      for (int64_t j = 0; j < dim; j++)
         this->m_basis[i][j] = m_basisOriginal[i][j];
}

//============================================================================

// Builds an upper-triangular primal basis in `dim` dimensions and computes its m-dual.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::buildDualBasis(int64_t dim) {
   buildBasis(dim);
   this->setDimDual(dim);
   mDualUpperTriangular(this->m_dualbasisOriginal, this->m_basisOriginal, this->m_modulo,
         this->m_maxDim);
   for (int64_t i = 0; i < dim; i++)
      for (int64_t j = 0; j < dim; j++)
         this->m_dualbasis[i][j] = m_dualbasisOriginal[i][j];
}

//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
// Same as in `MRGLatticeLac`.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::incDimBasis() {
   int64_t i, j, l;
   // assert(this->m_case1);  // This works only for case1.
   int64_t d = this->m_dim;
   if (m_dimM <= d) {
      m_M.SetDims(this->m_maxDim, this->m_maxDim);
      m_dimM = this->m_maxDim;
   }
   for (i = 0; i < d; i++) {
      for (j = 0; j < d; j++) {
         m_M[i][j] = this->m_basis[i][j];
         for (l = 0; l < j; l++) {
            m_M[i][j] -= m_M[i][l] * m_basisOriginal[l][j];
         }
         m_M[i][j] = m_M[i][j] / m_basisOriginal[j][j];
      }
   }
   for (i = 0; i < d; i++) {
      this->m_basis[i][d] = 0;
      for (j = 0; j < d; j++)
         MulAddTo(this->m_basis[i][d], m_M[i][j], m_basisOriginal[j][d]);
   }
   // Add last row from the stored primal basis.
   for (j = 0; j <= d; j++)
      this->m_basis[d][j] = m_basisOriginal[d][j];
   this->m_dim++;
   // std::cout << "m_basis after incDimBasis = \n" << this->m_basis << "\n";
}

//============================================================================

// Same as in 'MLCGLatticeLac'.
template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::incDimDualBasis() {
   int64_t d = this->m_dimdual;
   for (int i = 0; i < d; i++)
      this->m_dualbasis[i][d] = 0;
   for (int j = 0; j <= d; j++)
      this->m_dualbasis[d][j] = m_dualbasisOriginal[d][j];
   this->m_dimdual++;
}

template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   myExit("MLCGLatticeLac::buildProjection not yet implemented.");
}

template<typename Int, typename Real>
void MLCGLatticeLac<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   myExit("MLCGLatticeLac::buildProjectionDual not yet implemented.");
}

} // End namespace LatMRG
#endif


// This file is part of LatMRG.
//
// Copyright (C) 2012-2023  The LatMRG authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATMRG_MRGLATTICELAC_H
#define LATMRG_MRGLATTICELAC_H

#include <string>
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {

/**
 * This subclass of `MRGLattice` constructs and handles lattice bases built from MRGs as in `MRGLattice`,
 * but with arbitrary lacunary indices that can be spaced very far apart.
 * A special case of this is when they are regularly spaced by packets of the same size.
 *
 * To construct or increment the basis in that case, we proceed as described in Section 3.1.9 of the guide.
 * First, \f$P(z)\f$ must be set as the polynomial modulus in NTL.
 * Then for each lacunary index \f$\nu = i_j\f$, we compute the corresponding column by computing
 * \f$z^{\nu-1} \bmod P(z)\f$ using `power` from `polE` in NTL, then transforming its vector of coefficients
 * into the vector \f$(y_{\nu+k-2},\dots,y_{\nu-1})\f$ with the function `polyToColumn`.
 * This gives a set of generating vectors, which can be reduced to a basis, as in `MRGLattice`.
 * Under certain conditions, it is already a basis.
 */
template<typename Int, typename Real>
class MRGLatticeLac: public MRGLattice<Int, Real> {

public:

   /**
    * These four constructors have the same parameters as the parent constructors,
    * except for `maxCoord` which is not used.
    * We have the modulus \f$m\f$, vector of multipliers in `aa` (optional) or just the order \f$k\f$,
    * maximal dimension `maxDim` of a basis (primal or m-dual),
    * maximal dimension `maxDimProj` of a projection (for `buildProjection`), and the type of norm.
    * When not given, `maxDimProj` is taken equal to `maxDim`.
    * When `aa` is given, it must have `k+1` components with `a[j]`=\f$a_j\f$,
    * and `k` is deduced via `k = aa.length()-1`.
    * Otherwise, `aa` must be set afterwards via `setaa`.
    */
   MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, int64_t maxDimProj, NormType norm =
         L2NORM);

   MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   MRGLatticeLac(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj,
         NormType norm = L2NORM);

   MRGLatticeLac(const Int &m, int64_t k, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac(const MRGLatticeLac<Int, Real> &Lat);
   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac& operator=(const MRGLatticeLac<Int, Real> &Lat);
   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac() {
   }
   ;

   /**
    * Sets the vector of multipliers to `aa`.  This vector must have length \f$k+1\f$
    * where  \f$k\f$ is the order `k` of this MRG.
    * One will have \f$a_j=\f$`aa[j]` for \f$j=1,...,k\f$.
    * In contrast to the parent class, here no vector \f$\mathbf{y}\f$ is computed.
    */
   virtual void setaa(const IntVec &aa);

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
    * defined in the guide, for the current vector of multipliers and the current
    * set of lacunary indices \f$I = \{i_1,\dots,i_s\}\f$.
    * If some columns were already computed, it computes only the missing ones.
    * When the parameter `dim` is not given (with the default value of 0),
    * the number of columns is set to the number `m_s` of lacunary indices.
    * This function is called internally when needed, when we build a basis.
    */
   void computeVIks(int64_t dim = 0);

   /**
    * Builds an initial upper-triangular primal basis in `dim` dimensions.
    * This `dim` must be at least the order `k` and must not exceed `m_maxDim`.
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
    * This function uses the linear transformation method described in Section 3.1.9
    * of the guide.
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
    * Builds a basis for the projection of the current `MRGLatticeLac` onto the coordinates
    * in `proj` and puts it as the `m_basis` of `projLattice`.
    * This is done simply by projecting the current `m_basis`, as in the default
    * implementation in `IntLattice`.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` in `projLattice`.  Also done as the default in `IntLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;

   /**
    * Takes the polynomial `pcol` and returns in `col` the corresponding column in the
    * matrix of generating vectors.
    */
   void polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol);

protected:

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
    * The current number of columns for which `m_VIks` is computed.
    */
   int64_t m_dimVIks = 0;

   /**
    * This is the matrix \f$\mathbb{V}_{I,k\times s}\f$ defined in the guide.
    * It contains the ingredients we need to build a basis.
    * Each time the vector of coefficients \f$\mathbb{a}\f$ or the set \f$I\f$ of
    * lacunary indices changes, we must recompute this matrix.
    * This is done by `computeVIks`.
    */
   IntMat m_VIks;

   /**
    * Auxiliary matrix used to store a set of generating vectors used to compute a
    * triangular basis when `m_case == false`.
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

};

//===========================================================================
// Implementation:

template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj,
      NormType norm) :
      MRGLattice<Int, Real>(m, k, maxDim, norm) {
   this->m_order = k;
   this->m_aa.SetLength(k + 1);
   this->m_maxDimProj = maxDimProj;
   this->m_genTemp.SetDims(maxDimProj + k, maxDimProj);
   m_VIks.SetDims(k, maxDim);
   m_basisOriginal.SetDims(maxDim, maxDim);
   m_dualbasisOriginal.SetDims(maxDim, maxDim);
   m_M.SetDims(maxDim, maxDim);
   FlexModInt<Int>::mod_init(this->m_modulo);  // Sets `m` for polynomial arithmetic.
}

template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, int64_t k, int64_t maxDim, NormType norm) :
      MRGLatticeLac<Int, Real>(m, k, maxDim, maxDim, norm) {
}

template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      MRGLatticeLac<Int, Real>(m, aa.length() - 1, maxDim, maxDimProj, norm) {
   this->setaa(aa);
}

template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim,
      NormType norm) :
      MRGLatticeLac<Int, Real>(m, aa, maxDim, maxDim, norm) {
}

/*
 //============================================================================

 template<typename Int, typename Real>
 MRGLatticeLac<Int,Real>::MRGLatticeLac(const MRGLatticeLac<Int, Real> &Lat) {
 this->m_modulo = Lat.m_modulo;
 this->m_maxDim = Lat.m_maxDim;
 this->m_dim = Lat.m_dim;
 this->m_dimdual = Lat.m_dimdual;
 this->m_norm = Lat.m_norm;
 this->m_vecNorm = RealVec(Lat.m_vecNorm);
 this->m_dualvecNorm = RealVec(Lat.m_dualvecNorm);
 this->m_basis = Lat.m_basis;
 this->m_dualbasis = Lat.m_dualbasis;
 setLac(Lat.m_lac, false);
 FlexModInt<Int>::mod_init(this->m_modulo);
 setaa(Lat.m_aa, true);
 }

 //============================================================================

 // Copy Constructor
 template<typename Int, typename Real>
 MRGLatticeLac<Int, Real> & MRGLatticeLac<Int, Real>::operator= (const MRGLatticeLac<Int, Real> & Lat) {
 this->m_modulo = Lat.m_modulo;
 this->m_maxDim = Lat.m_maxDim;
 this->m_dim = Lat.m_dim;
 this->m_dimdual = Lat.m_dimdual;
 this->m_norm = Lat.m_norm;
 this->m_vecNorm = RealVec(Lat.m_vecNorm);
 this->m_dualvecNorm = RealVec(Lat.m_dualvecNorm);
 this->m_basis = Lat.m_basis;
 this->m_dualbasis = Lat.m_dualbasis;
 setLac(Lat.m_lac, false);
 FlexModInt<Int>::mod_init(this->m_modulo);
 setaa(Lat.m_aa, true);
 return *this;
 }
 */

//============================================================================
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::setLac(const IntVec &lac) {
   int64_t k = this->m_order;
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
}

//============================================================================

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::setaa(const IntVec &aa) {
   assert(aa.length() == this->m_order + 1);  // `aa` must have the right length.
   this->m_aa = aa;
   m_dimVIks = 0;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   // Set the characteristic polynomial of the recurrence in NTL.
   typename FlexModInt<Int>::PolX m_Pz;
   for (int j = 1; j < this->m_aa.length(); j++) {
      SetCoeff(m_Pz, this->m_aa.length() - j - 1, FlexModInt<Int>::to_Int_p(-this->m_aa[j]));
   }
   SetCoeff(m_Pz, this->m_aa.length() - 1, 1);
   FlexModInt<Int>::PolE::init(m_Pz);
}

//============================================================================

// This function computes the polynomials z^{\mu-1} mod P(z) for \mu in the set of
// lacunary indices and use them to compute the first `dim` columns of the matrix `m_VIks`,
// as in Section 3.1.9 of the guide.
// If some columns were already computed, it computes only the missing ones.
// By default, `dim` is set to the number of lacunary indices.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::computeVIks(int64_t dim) {
   assert(dim <= min(this->m_maxDim, m_s));  // Upper bound on `dim`.
   if (dim == 0) dim = m_s;
   if (dim <= m_dimVIks) return;     // Nothing to do.
   int64_t k = this->m_order;
   int64_t i, j;
   IntVec col;
   typename FlexModInt<Int>::PolE polz, polPower;  // polynomials z and z^{\nu-1} mod P(z).
   std::string str = "[0 1]";
   std::istringstream in(str);
   in >> polz;

   // Fill the rows of VIks, column by column.
   for (j = m_dimVIks; j < dim; j++) {
      if (m_case1 && j < k) for (i = 0; i < k; i++)
         m_VIks[i][j] = (i == j);
      // If index j has just increased by 1, the new column is easy to compute.
      else if ((j > 0) && (m_lac[j] == m_lac[j - 1] + 1)) {
         for (i = 1; i < k; i++)
            m_VIks[i][j] = m_VIks[i - 1][j - 1];
         m_VIks[0][j] = 0;
         for (int64_t jj = 1; jj <= k; jj++)
            m_VIks[0][j] += this->m_aa[jj] * m_VIks[jj - 1][j - 1] % this->m_modulo;
      }
      // Otherwise compute the polynomial power.
      else {
         power(polPower, polz, m_lac[j] - 1);
         this->polyToColumn(col, polPower);
         for (i = 0; i < k; i++) {
            m_VIks[i][j] = col[i];
         }
      }
   }
   m_dimVIks = dim;
}

//============================================================================

// Builds an upper-triangular `basis` in `dim` dimensions, using `m_VIks`.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildBasisOriginal(IntMat &basis, int64_t dim) {
   int64_t k = this->m_order;
   assert(dim >= k);
   assert(dim <= min(this->m_maxDim, m_s));
   // this->setDim(dim);
   if (dim > m_dimVIks) computeVIks(dim);
   int64_t i, j;
   // Fill the first k rows using m_VIks.
   IntMat & pbasis = basis;      // reference to a basis.
   if (!m_case1) pbasis = this->m_genTemp;  // Should be rare.
   for (i = 0; i < k; i++)
      for (j = 0; j < dim; j++)
         pbasis[i][j] = m_VIks[i][j];
   // Fill the other rows.
   if (m_case1) for (i = k; i < dim; i++)
      for (j = 0; j < dim; j++)
         pbasis[i][j] = (i == j) * this->m_modulo;
   else {
      for (i = 0; i < dim; i++)
         for (j = 0; j < dim; j++)
            pbasis[i + k][j] = (i == j) * this->m_modulo;
      upperTriangularBasis(this->m_basis, pbasis, this->m_modulo, dim + k, dim);
   }
}

//============================================================================

// Builds an upper-triangular `basis` in `dim` dimensions, using `m_VIks`.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildBasis(int64_t dim) {
   this->m_dim = dim;
   if (m_dimVIks < this->m_maxDim) computeVIks(this->m_maxDim);
   buildBasisOriginal(m_basisOriginal, this->m_maxDim);
   for (int64_t i = 0; i < dim; i++)
      for (int64_t j = 0; j < dim; j++)
         this->m_basis[i][j] = m_basisOriginal[i][j];
}

//============================================================================

// Builds an upper-triangular primal basis in `dim` dimensions and computes its m-dual.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildDualBasis(int64_t dim) {
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
// The algorithm described in Section 3.1.9 of the guide is implemented.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimBasis() {
   int64_t i, j, l;
   assert(this->m_case1);  // This works only for case1.
   this->m_dim++;
   int64_t d = this->m_dim;

   // IntMat M;
   // Calculate the matrix M according to guide. It stores the transformation from the standard basis to the current one.
   // M.SetDims(d-1, d-1);  // This allocates space for a new matrix object from scratch each time we call this function!
   for (i = 0; i < d - 1; i++) {
      for (j = 0; j < d - 1; j++) {
         m_M[i][j] = this->m_basis[i][j];
         for (l = 0; l < j; l++) {
            m_M[i][j] -= m_M[i][l] * m_basisOriginal[l][j];
         }
         m_M[i][j] = m_M[i][j] / m_basisOriginal[j][j];
      }
   }
   // Calculate the new last column by applying M to the last column of the stored primal basis.
   //IntMat copy_curr_column, new_last_column;
   //copy_curr_column.SetDims(d-1,1);  // Again another matrix creation!

   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = 0;
      for (j = 0; j < d - 1; j++)
         MulAddTo(this->m_basis[i][d - 1], m_M[i][j], m_basisOriginal[j][d - 1]);
   }
   // Add last row from the stored primal basis.
   for (j = 0; j < d; j++)
      this->m_basis[d - 1][j] = m_basisOriginal[d - 1][j];
}

//============================================================================

// Very similar to the implementation in 'MRGLattice' except that the copy
// of the dual basis has already been built upon creation of the object.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimDualBasis() {
   this->m_dimdual++;
   int64_t d = this->m_dimdual;
   // Put zeros in new column and add a new row.
   for (int i = 0; i < d - 1; i++)
      this->m_dualbasis[i][d - 1] = 0;
   for (int j = 0; j < d; j++)
      this->m_dualbasis[d - 1][j] = m_dualbasisOriginal[d - 1][j];
}

//============================================================================

/*
 template<typename Int, typename Real>
 bool MRGLatticeLac<Int, Real>::buildProjection(IntMat &basis, int64_t dimbasis, IntMat &pbasis,
 const Coordinates &proj) {
 projectionConstructionUpperTri<Int>(pbasis, basis, proj, this->m_modulo,
 this->m_dim);
 return false;
 }
 */

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   myExit("MRGLatticeLac::buildProjection not yet implemented.");
}

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   myExit("MRGLatticeLac::buildProjectionDual not yet implemented.");
}

//============================================================================

// This function applies phi inverse as described in Section 3.1.2 of the guide, and reverses the coordinates.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol) {
   int i, j, k;
   k = this->m_order;
   IntVec c;
   c.SetLength(k);   // Allocates memory for a new vector!
   col.SetLength(k);
   for (j = 1; j < k + 1; j++) {
      c[j - 1] = 0;
      NTL::conv(c[j - 1], coeff(rep(pcol), k - j));
      for (i = 1; i < j; i++) {
         c[j - 1] += this->m_aa[i] * c[j - 1 - i];
      }
      c[j - 1] = c[j - 1] % this->m_modulo;
   }
   for (j = 0; j < k; j++) {
      col[j] = c[k - j - 1];
   }
}

}
#endif

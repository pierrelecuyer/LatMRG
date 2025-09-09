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

#ifndef LATMRG_LCGLATTICELAC_H
#define LATMRG_LCGLATTICELAC_H

#include <string>
#include "NTL/ZZ.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/LCGLattice.h"

namespace LatMRG {

/**
 * This subclass of `LCGLattice` constructs and handles lattice bases built from LCGs as in `LCGLattice`,
 * but with arbitrary lacunary indices that can be spaced very far apart.
 * A special case of this is when they are regularly spaced by packets of the same size.
 * No need for polynomial arithmetic in this class, just exponentiation mod m.
 */
template<typename Int, typename Real>
class LCGLatticeLac: public LCGLattice<Int, Real> {

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
   LCGLatticeLac(const Int &m, const Int &a, int64_t maxDim, int64_t maxDimProj, NormType norm =
         L2NORM);

   LCGLatticeLac(const Int &m, const Int &a, int64_t maxDim, NormType norm = L2NORM);

   LCGLatticeLac(const Int &m, int64_t maxDim, int64_t maxDimProj,
         NormType norm = L2NORM);

   LCGLatticeLac(const Int &m, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   // LCGLatticeLac(const LCGLatticeLac<Int, Real> &Lat);
   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   // LCGLatticeLac& operator=(const LCGLatticeLac<Int, Real> &Lat);
   /**
    * Destructor.
    */
   virtual ~LCGLatticeLac() {
   }
   ;

   /**
    * Sets the multiplier to `a`.
    */
   virtual void seta(const Int &a);

   /**
    * Sets the vector of lacunary indices to `lac`, whose length should
    * not exceed `m_maxDim`. This length is put into `m_s`.
    * The indices are assumed to be in increasing order.
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
    * Compute the values \f$(y_0, y_{i_1}, ..., y_{i_{t-1}})\f$ for \f$\t=\f$`dim`,
    * which cannot exceed `m_maxDim` or `m_s`.
    * If `dim = 0`, the function takes `dim` equal to `m_s`, the number of lacunary indices.
    * If some values are already computed, it computes only the missing ones, if any.
    */
   virtual void computeyLac(int64_t dim = 0);

   /**
    * Builds an initial upper-triangular primal basis in `dim` dimensions.
    * This `dim` must not exceed `m_maxDim`.
    * The function `computeyLac` is called internally if needed.   ????
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
    * Builds a basis for the projection of the current `LCGLatticeLac` onto the coordinates
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


protected:

   /**
    * Stores the set \f$I\f$ of lacunary indices and its size \f$s\f$.
    */
   IntVec m_lac;

   int64_t m_s = 0;

   /**
    * This vector corresponds to the first row of the matrix \f$\mathbb{V}_{I,k\times s}\f$
    * defined in the guide, which is `m_VIks` in `MRGLatticeLac`.
    * It contains the ingredients we need to build a basis.
    * Each time the vector of coefficients \f$\mathbb{a}\f$ or the set of
    * lacunary indices changes, we must recompute this vector.
    * This is done by `computeyLac`.
    */
   IntVec m_yLac;

   /**
    * The current number of columns for which `m_y` is computed.
    */
   int64_t m_dimyLac = 0;

   IntMat m_M;  // Matrix M used in `incDimBsis`.
   int64_t m_dimM = 0;  // Its current dimension.
};

//===========================================================================
// Implementation:

template<typename Int, typename Real>
LCGLatticeLac<Int, Real>::LCGLatticeLac(const Int &m, int64_t maxDim, int64_t maxDimProj,
      NormType norm) :
      LCGLattice<Int, Real>(m, maxDim, norm) {
   this->m_maxDimProj = maxDimProj;
   m_yLac.SetLength(maxDim);
   m_M.SetDims(maxDim, maxDim);
}

template<typename Int, typename Real>
LCGLatticeLac<Int, Real>::LCGLatticeLac(const Int &m, int64_t maxDim, NormType norm) :
      LCGLatticeLac<Int, Real>(m, maxDim, maxDim, norm) {
}

template<typename Int, typename Real>
LCGLatticeLac<Int, Real>::LCGLatticeLac(const Int &m, const Int &a, int64_t maxDim,
      int64_t maxDimProj, NormType norm) :
      LCGLatticeLac<Int, Real>(m, maxDim, maxDimProj, norm) {
   this->seta(a);
}

template<typename Int, typename Real>
LCGLatticeLac<Int, Real>::LCGLatticeLac(const Int &m, const Int &a, int64_t maxDim,
      NormType norm) :
      LCGLatticeLac<Int, Real>(m, a, maxDim, maxDim, norm) {
}

//============================================================================
template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::setLac(const IntVec &lac) {
   // std::cout << "\n setLac, k = " << k << ",  maxDim = " << this->m_maxDim << ",\n lac = " << lac << "\n";
   assert((lac.length() >= 1) && (lac.length() <= this->m_maxDim));
   m_lac = lac;
   m_s = lac.length();
   m_dimyLac = 0;    // The vector y must be recomputed.
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
}

//============================================================================

template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::seta(const Int &a) {
   this->m_a = a;
   m_dimyLac = 0;    // The vector y must be recomputed.
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
}

//============================================================================

// Computes the powers of `a` modulo `m` to fill `m_yLac` up to `dim` dimensions.
template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::computeyLac(int64_t dim) {
   assert(dim <= m_s);   // Upper bound on `dim`.
   if (dim == 0) dim = m_s;
   if (dim <= m_dimyLac) return;     // Nothing to do.
   int64_t j;
   m_yLac[0] = 1;
   for (j = max(m_dimyLac, 1); j < dim; j++) {
      // if (j == 0) m_yLac[j] = 1;
      if ((m_lac[j] == m_lac[j - 1] + 1))
         m_yLac[j] = this->m_a * m_yLac[j - 1] % this->m_modulo;
      else
         PowerMod(m_yLac[j], this->m_a, m_lac[j]-1, this->m_modulo);
   }
   m_dimyLac = dim;
}

//============================================================================

// Builds an upper-triangular `basis` in `dim` dimensions, using `m_VIks`.
template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::buildBasis(int64_t dim) {
   this->setDim(dim);
   int64_t i, j;
   if (dim > m_dimyLac) computeyLac(); // Compute for m_s coordinates.
   for (j = 0; j < dim; j++)
      this->m_basis[0][j] = m_yLac[j];
   for (i = 1; i < dim; i++)
      for (j = 0; j < dim; j++)
         this->m_basis[i][j]  = (i == j) * this->m_modulo;
}

//============================================================================

template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::buildDualBasis(int64_t dim) {
   this->setDimDual(dim);
   if (dim > m_dimyLac) computeyLac();
   this->m_dualbasis[0][0] = this->m_modulo;
   int64_t i, j;
   for (i = 1; i < dim; i++)
      this->m_dualbasis[i][0] = -this->m_yLac[i] % this->m_modulo;
   for (i = 0; i < dim; i++)
      for (j = 1; j < dim; j++)
         this->m_dualbasis[i][j] = (i == j);
}

//============================================================================

template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::incDimBasis() {
   int64_t i, j;
   int64_t d = this->m_dim;
   if (d >= m_dimyLac) computeyLac(d+1);
   // Here we reserve space for M, only if it is needed.
   if (m_dimM <= d) {
      m_M.SetDims(this->m_maxDim, this->m_maxDim);
      m_dimM = this->m_maxDim;
   }
   for (i = 0; i < d; i++) {
      m_M[i][0] = this->m_basis[i][0];
      for (j = 1; j < d; j++)
         m_M[i][j] = (this->m_basis[i][j] - m_M[i][0] * m_yLac[j-1]) / this->m_modulo;
   }
   // To obtain the new last column, multiply by m_M the original last column
   // without its last coordinate, which is (m_yLac[d], 0, ..., 0).
   for (i = 0; i < d; i++) {
      this->m_basis[i][d] = m_M[i][0] * m_yLac[d];
   }
   this->m_dim++;
}

//============================================================================

// Add one zero coordinate to each vector, and add the last row of W_I^{(p)}.
template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::incDimDualBasis() {
   int64_t d = this->m_dimdual;
   if (d >= m_dimyLac) computeyLac(d+1);
   // Put zeros in new column and add a new row.
   for (int i = 0; i < d; i++)
      this->m_dualbasis[i][d] = 0;
   this->m_dualbasis[d][0] = -m_yLac[d];
   this->m_dualbasis[d][d] = 1;
   for (int j = 1; j < d; j++)
      this->m_dualbasis[d][j] = 0;
   this->m_dimdual++;
}

//============================================================================

template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   myExit("LCGLatticeLac::buildProjection not yet implemented.");
}

template<typename Int, typename Real>
void LCGLatticeLac<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   myExit("LCGLatticeLac::buildProjectionDual not yet implemented.");
}

}
#endif

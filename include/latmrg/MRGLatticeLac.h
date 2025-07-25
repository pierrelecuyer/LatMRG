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
   // Parent:
   // MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Constructor with modulus of congruence \f$m\f$, order of the recurrence
    * \f$k\f$, multipliers in `aa`, and maximal dimension `maxDim`.
    * The length of basis vectors is computed with `norm`.
    * The basis is built for the lacunary indices in `lac`.
    * The vector `aa` must have k+1 components with `a[j]`=\f$a_j\f$.
    */
   MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, IntVec &lac,
         NormType norm = L2NORM);
   
   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac(const MRGLatticeLac<Int, Real> &Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac& operator=(const MRGLatticeLac<Int, Real> &Lat);

   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac(){ };

   /**
    * Sets the lacunary indices for this lattice to `lac`.
    * If 'buildBasisCopy' is set to true, then 'm_copy_primal_basis' and
    * 'm_dual_copy_basis' are updated according to the new lacunary indices.
    */
   void setLac(const IntVec &lac, bool buildBasisCopy = true);

   
   /**
    * Sets the generating vector to `aa`.
    * If 'buildBasisCopy' is set to true, then 'm_copy_primal_basis' and
    * 'm_dual_copy_basis' are updated according to the new generator vector.
    */
   void setaa(const IntVec &lac, bool buildBasisCopy = true) override;

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   Int& getLac(int j);

protected:


   /**
    * This function overrides the correpsonding protected function in 'MRGLattice'.
    * It Builds a basis directly in `d` dimensions, as explained in Section 3.1.9 of
    * the LatMRG guide.  Must have d <= m_maxdim. The basis matrix is taken as a 
    * parameter.
    */ 
   void buildBasis0Pol(IntMat &basis, int64_t d) override;
   
   /**
    * This functions overrides the corresponding protected function in 'MRGLattice',
    * which is necessary because the entire copy of the primal basis is build upon
    * creation of the object.
    */
   void buildDualBasis0(IntMat &basis, int64_t d) override;
   
   /**
    * This function overrides the correpsonding protected function in 'MRGLattice'.
    * It increases the dimension of the given basis from d-1 to d dimensions.
    * One new column is calclated using the algorithm from the guide.
    */
   void incDimBasis0(IntMat &basis, int64_t d) override;

   /**
    * This function overrides the correpsonding protected function in 'MRGLattice'.
    * It increases the dimension of a given dual basis from d-1 to d dimensions.
    * One new column is calclated using the polynomial representation.
    */
   void incDimDualBasis0Pol(IntMat &basis, int64_t d) override;
      
   /**
    * Stores the lacunary indices.
    */
   IntVec m_lac;

};


//===========================================================================
// Implementation:


// Constructor.
template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim,
      IntVec &lac, NormType norm) : MRGLattice<Int, Real>(m, aa, maxDim, norm) {
      this->m_modulo = m;
      this->m_maxDim = maxDim;   
      this->m_dim = 0;
      this->m_copy_primal_basis.SetDims(maxDim, maxDim);
      this->m_copy_dual_basis.SetDims(maxDim, maxDim);
      setLac(lac, false);
      this->setaa(aa, false);
      FlexModInt<Int>::mod_init(m);
      this->buildyPol(maxDim + this->m_order - 1);
      // Immediately build a copy of the full basis and store in copy_primal
      buildBasis0Pol(this->m_copy_primal_basis, maxDim);
      mDualUpperTriangular(this->m_copy_dual_basis, this->m_copy_primal_basis, this->m_modulo, maxDim);
      
}


//============================================================================

// Copy Constructor
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
   setaa(Lat.m_aCoeff, true);
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
   setaa(Lat.m_aCoeff, true);
   return *this;
}

//============================================================================

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::setLac(const IntVec &lac, bool buildBasisCopy) {
   m_lac = lac;
   if (buildBasisCopy) {
      buildBasis0Pol(this->m_copy_primal_basis, this->m_maxDim);
      mDualUpperTriangular(this->m_copy_dual_basis, this->m_copy_primal_basis, this->m_modulo, this->m_maxDim);
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::setaa(const IntVec &aa, bool buildBasisCopy) {
   this->m_aCoeff = aa;
   this->m_order = aa.length() - 1;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   if (buildBasisCopy) {
      buildBasis0Pol(this->m_copy_primal_basis, this->m_maxDim);
      mDualUpperTriangular(this->m_copy_dual_basis, this->m_copy_primal_basis, this->m_modulo, this->m_maxDim);
   }
}

//============================================================================

// Replaces the corresponding function in `MRGLattice`.

// Builds a basis directly in `d` dimensions, as explained in Section 3.1.9 of
// the guide of LatMRG.  Must have d <= m_maxdim.
// The implementation is very except that the indices are usually far apart.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildBasis0Pol(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   assert(d <= m_lac.length());
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   int64_t i, j;
   IntVec col;
   
   typename FlexModInt<Int>::PolE polDegOne;
   typename FlexModInt<Int>::PolE polPower;
      
   // Auxilliary variables to calculate powers p^n(z)
   std::string str = "[0 1]";
   std::istringstream in (str);
   in >> polDegOne;   
   // Calculate powers p^\mu-1 for \mu in the set of lacunary indices
   // And fill the first k rows
   for (j = 0; j < d; j++) {
      power(polPower, polDegOne, getLac(j) - 1);
      this->polyToColumn(col, polPower);
      for (i = 0; i < dk; i++) {
         basis[i][j] = col[i];
      }
   }
   // Fill the rest of the rows
   for (i = dk; i < d; i++) { 
      for (j = 0; j < d; j++)
         basis[i][j] = (i == j) * this->m_modulo;
   }
}

//============================================================================

// The implemtnation needs to take into account that 'm_copy_primal_basis' is built
// upon creation of the object.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildDualBasis0(IntMat &basis, int64_t d) {
      mDualUpperTriangular(basis, this->m_copy_primal_basis, this->m_modulo, d);
}


//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
// The algorithm described in Section 3.1.9 of the guide is implemented.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   int64_t i, j, l;
   IntMat M; 
   // Calculate the matrix M according to guide. It stores the transformation from the standard basis to the current one.
   M.SetDims(d-1, d-1);
   for (i = 0; i < d-1; i++) {
      for (j = 0; j < d-1; j++) {
         M[i][j] = basis[i][j];
            for (l = 0; l < j; l++) {
               M[i][j] -= M[i][l]*this->m_copy_primal_basis[l][j];
            }
         M[i][j] = M[i][j] / this->m_copy_primal_basis[j][j];
      }
   }
   
   // Calculate the new last column by applying M to the last column of the stored primal basis.
   IntMat copy_curr_column, new_last_column;
   copy_curr_column.SetDims(d-1,1);
   for (i = 0; i < d-1; i++)  
      copy_curr_column[i][0] = this->m_copy_primal_basis[i][d-1];
   mul(new_last_column, M, copy_curr_column);
   for (i = 0; i < d-1; i++)
      basis[i][d-1] = new_last_column[i][0];
   
   // Add last row from the stored primal basis.
   for (j = 0; j < d; j++)
      basis[d-1][j] = this->m_copy_primal_basis[d-1][j];
}

//============================================================================

// Very similar to the implementation in 'MRGLattice' except that the copy
// of the dual basis has already been built upon creation of the object.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimDualBasis0Pol(IntMat &basis, int64_t d) {
   // Add zeros to last column
   for (int i = 0; i < d-1; i++)
      basis[i][d-1] = 0;
   // Add last row
   for (int i = 0; i < d; i++)
      basis[d-1][i] = this->m_copy_dual_basis[d-1][i];
}

}
#endif

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
//#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {


/**
 * This subclass of `MRGLattice` constructs and handles lattice bases built from MRGs as in `MRGLattice`,
 * but with arbitrary lacunary indices that can be spaced very far apart.
 * A special case of this is when they are regularly spaced by packets of the same size.
 *
 * To construct or increment the basis in that case, we proceed as described in Section 4.1.9 of the guide.
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
    */
   void setLac(const IntVec &lac) {m_lac = lac;};

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   Int& getLac(int j);

protected:

   /**
    * This function overrides the correpsonding protected function in 'MRGLattice'.
    * It Builds a basis directly in `d` dimensions, as explained in Section 4.1.9 of
    * the LatMRG guide.  Must have d <= m_maxdim. The basis matrix is taken as a 
    * parameter.
    */ 
   void buildBasis0Pol(IntMat &basis, int64_t d) override;
   
   /**
    * This function overrides the correpsonding protected function in 'MRGLattice'.
    * It increases the dimension of given basis from d-1 to d dimensions.
    * One new column is calclated using the polynomial representation.
    */
   void incDimBasis0(IntMat &basis, int64_t d) override;
   
   /**
    * The lacunary indices.
    */
   IntVec m_lac;

};


//===========================================================================
// Implementation:


// Constructor.
template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim,
      IntVec &lac, NormType norm) : MRGLattice<Int, Real>(m, aa, maxDim, norm) {
      setLac(lac);
      this->setaa(aa);
      FlexModInt<Int>::mod_init(m);
      this->buildyPol(maxDim + this->m_order - 1);
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
   setLac(Lat.m_lac);
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
   setLac(Lat.m_lac);
   return *this;
}

//============================================================================

// Replaces the corresponding function in `MRGLattice`.

// Builds a basis directly in `d` dimensions, as explained in Section 4.1.9 of
// the guide of LatMRG.  Must have d <= m_maxdim.
// This should be very similar to `buildProjection0`, except that the indices are usually far apart.

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
      power(polPower, polDegOne, m_lac[j] - 1);
      this->polyToColumn(col, polPower);
      for (i = 0; i < dk; i++)
         basis[i][j] = col[i];
   }
   // Fill the rest of the rows
   for (i = dk; i < d; i++) { 
      for (j = 0; j < d; j++)
         basis[i][j] = (i == j) * this->m_modulo;
   }
}


//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
// We compute one new column using the polynomial representation.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   IntVec col;
   int64_t i;
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   typename FlexModInt<Int>::PolE polDegOne;
   typename FlexModInt<Int>::PolE polPower;
   assert(d <= this->m_maxDim);
   // Auxilliary variables to calculate powers p^n(z)
   std::string str = "[0 1]";
   std::istringstream in (str);
   in >> polDegOne;   
   power(polPower, polDegOne, m_lac[d-1] - 1);
   this->polyToColumn(col, polPower);
   // Put the coefficients in the desired row
   for (i = 0; i < dk; i++)
      basis[i][d-1] = col[i];
   // Fill the rest of the row
   for (i = dk; i < d; i++)
      basis[i][d-1] = (i == d-1) * this->m_modulo;
}

}
#endif

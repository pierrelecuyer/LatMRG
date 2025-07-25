// This file is part of LatMRG.
//
// Copyright (C) 2012-2022  The LatMRG authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit? de Montr?al.
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

#ifndef LATMRG_MRGLATTICE_H
#define LATMRG_MRGLATTICE_H

#include <cassert>
#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FlexTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/BasisConstruction.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {

/**
 * This subclass of `IntLatticeExt` defines an MRG lattice and is similar to Rank1Lattice,
 * but constructs lattices associated with multiple recursive generators (MRGs) with 
 * modulus m, order k, and vector of multipliers a = (a_1, . . . , a_k).
 *
 * The functions `buildBasis` and `incDimBasis` always build and update the primal basis.
 * To build the m-dual, use `buildDualBasis`.
 */

template<typename Int, typename Real>
class MRGLattice: public IntLatticeExt<Int, Real> {

private:
   // typedef NTL::vector<Int> IntVec;
   // typedef NTL::matrix<Int> IntMat;
   // typedef NTL::vector<Real> RealVec;

public:

   /**
    * This constructor takes as input the modulus `m`, the vector of multipliers `aa`,
    * and the norm used to measure the vector lengths.
    * The vector `aa` must have length `k+1`, with \f$a_j\f$ in `aa[j]`.
    * The maximal dimension `maxDim` will be the maximal dimension of the basis.
    * This constructor does not build the basis, so we can build it for a smaller
    * number of dimensions or only for selected projections.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /*
    * Copy constructor.
    */
   MRGLattice(const MRGLattice<Int, Real> &Lat);

   /**
    * Assigns `Lat` to this object.
    */
   MRGLattice& operator=(const MRGLattice<Int, Real> &Lat);

   /**
    * Destructor.
    */
   ~MRGLattice();

   /**
    * Sets the vector of multipliers. The order `k` of the lattice is set equal to
    * the length of this vector, minus 1. `aa[j]` must contain \f$a_j\f$ for j=1,...,k.
    */
   virtual void setaa(const IntVec &aa);   
  
   /**
    * Builds a basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
    * This initial primal basis will be upper triangular -  based on V^{(0)} or V^{(p)}
    * depending on which approach is chosen by the user.
    */
   void buildBasis(int64_t dim);

   /**
    * Builds the m-dual lower triangular basis directly in `dim` dimensions. This `dim` must 
    * not exceed `maxDim`. 
    * In order to build the m-dual basis, the primal basis is needed, even if it has not been
    * built yet. Therefore, a copy of the basis matrix is always created upon building the 
    * m-dual basis. It is stored in the variable 'm_copy_primal_basis'.
    */
   void buildDualBasis(int64_t dim); 

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `maxDim`.
    */
   void incDimBasis();

   /**
    * Increases the current dimension of the m-dual basis by 1.
    * The new increased dimension must not exceed `maxDim`.
    * This function uses the simplified method for MRG lattices given in the lattice tester guide.
    * It requires to increase the dimension of 'm_copy_primal_basis' as well.
    */
   void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. A basis for the projection of this
    * `MRGLattice` over the coordinates in `proj` is returned in `projLattice`.
    * The implementation used here exploits the rank-k lattice structure and it
    * is faster than the general one. See Section 3.1.7 of the LatMRG guide.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Overrides the same function from `IntLattice`.
    * When all the first k coordinates {1,...,k} belong to the projection, both the primal
    * and m-dual constructions are direct, just by selecting the rows and columns
    * whose indices are in `proj`. Otherwise an upper triangular basis is constructed
    * for the basis and a lower-triangular basis for the m-dual.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /** 
    * If 'usepol' is set to true then the polynomial basis V^{(p)}, see Section 3.1.5 of the guide,
    * is used. Otherwise the standard basis V^{(0)} is used, see Section 3.1.4 of the guide.
    */
   void setUsePolynomialBasis (const bool usepol) { use_polynomial_basis = usepol; }

   /**
    * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
    * where `dim` is the current lattice dimension.
    */
   std::string toStringCoef() const;
   
protected:
  
   /**
    * The following protected functions take the basis as a parameter for more flexibility.
    * They are used inside buildBasis, buildBasisDual, incDimBasis, etc.
    * For building the basis either the standard approach or the polynomial approach can be chosen.
    * The polynomial approach is implemented in those functions ending with Pol
    */
   virtual void buildBasis0(IntMat &basis, int64_t d);
   
   virtual void buildBasis0Pol(IntMat &basis, int64_t d);
   
   virtual void buildDualBasis0(IntMat &basis, int64_t d);      

   virtual void incDimBasis0(IntMat &basis, int64_t d);
   
   virtual void incDimDualBasis0(IntMat &basis, int64_t d);   
   
   virtual void incDimDualBasis0Pol(IntMat &basis, int64_t d);

   bool buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis, const Coordinates &proj);

   /**
    * Takes the polynomial `pcol` and returns in `col` the corresponding column in the
    * matrix of generating vectors.
    */
   void polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol);

   /**
    * Builds the vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used to build the matrix V^{(p)},
    * see Section 3.1.2 of the guide for details.
    */
   void buildyPol(int64_t dim);   

   /**
    * A vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used in the matrix V^{(p)}.
    */
   IntVec m_y;

   /**
    * This auxiliary matrix is used to store the generating vectors of a projections
    * before reducing them into a triangular basis.
    */
   IntMat m_genTemp;
   
   /**
    * For generating the dual basis or increasing its dimension, we need a copy of the
    * the primal basis.
    */
   IntMat m_copy_primal_basis;
   
   /**
    * If we want to increase the dimension of the dual basis with the polynomial approach
    * we also need a copy of the dual basis.
    */   
   IntMat m_copy_dual_basis;

   /**
    * Boolean variable which decides whether the polynomial basis V^{(p)} is used instead of V^{(0)}
    */   
   bool use_polynomial_basis = true;

   // Order of this MRG.
   int m_order;

   /**
    * The coefficients \f$a_1, ..., a_k\f$ of the MRG recurrence, a_j stored in `m_aCoeff[j]`.
    */
   IntVec m_aCoeff;


};

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_modulo = m;
   this->m_maxDim = maxDim;   
   setaa(aa);
   this->m_dim = 0;
   m_genTemp.SetDims(maxDim, maxDim); 
   m_copy_primal_basis.SetDims(maxDim, maxDim);
   m_copy_dual_basis.SetDims(maxDim, maxDim);   
   FlexModInt<Int>::mod_init(m);
   buildyPol(maxDim + m_order - 1);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::~MRGLattice() {
   m_aCoeff.kill();
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(const MRGLattice<Int, Real> &lat) {
   if (this == &lat) return *this;
   // this->copy(lat); CW: copy constructor currently not available
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_y = lat.m_y;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.getNormType()) {
   // this->copy(lat); CW: copy constructor currently not available
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_y = lat.m_y;
}

//============================================================================

/**
 * Sets the generating vector to `aa`.
 */
template<typename Int, typename Real>
void MRGLattice<Int, Real>::setaa(const IntVec &aa) {
   m_aCoeff = aa;
   m_order = aa.length() - 1;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildyPol(int64_t dim) {
   int64_t k = m_order;
   int64_t j, i;
   int64_t n;
   typename FlexModInt<Int>::PolE polDegOne;
   typename FlexModInt<Int>::PolE polPower;     
   typename FlexModInt<Int>::PolX m_Pz;  
   IntVec col;
   
   // Set the characteristic polynomial of the recurrence   
   for (j = 1; j < this->m_aCoeff.length(); j++) {
     SetCoeff(m_Pz, this->m_aCoeff.length() - j - 1, FlexModInt<Int>::to_Int_p(-this->m_aCoeff[j]));
   }   
   SetCoeff(m_Pz, this->m_aCoeff.length() - 1, 1);    
   FlexModInt<Int>::PolE::init(m_Pz);   
   
   // Auxilliary variables to calculate powers p^n(z)
   std::string str = "[0 1]";
   std::istringstream in (str);
   in >> polDegOne;   
  
   // Builde the vector y according to the algorithm described in Section 3.1.2 of the guide.
   // The vector is stored in the variable 'm_y'.
   m_y.SetLength(dim);
   for (j = 0; j < dim; j++)
     m_y[j] = 0;
   m_y[k-1] = 1;
   n = ceil(dim/k);
   for (j = 1; j < n+1; j++) {
      // Calculate powers p^\mu-1 
      power(polPower, polDegOne, j*k);
      polyToColumn(col, polPower);
      // And fill the next k entries
      for (i = 0; i < k; i++) {
         if (j*k+i < dim) {
            m_y[j*k+i] = col[k-i-1];
         }
      }
   }   
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension isset to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   this->setDim(d);
   if (use_polynomial_basis)
      this->buildBasis0Pol(this->m_basis, d);
   else
      this->buildBasis0(this->m_basis, d);
   this->setNegativeNorm();
}

//============================================================================

// Builds the upper-triangular basis V^{(0)} directly in `d` dimensions, as explained in Section 3.1.4 of
// the guide of LatMRG, puts this matrix in `basis`.  Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = m_order;
   int64_t dk = min(d, k);
   int64_t i, j, jj;
   // Put the lower-triangular part of the identity matrix in the upper left corner
   for (i = 0; i < dk; i++) {
      for (j = 0; j <= i; j++)
         basis[i][j] = (i == j);  // Avoid "if" statements.
   }
   // Put m times the identity matrix to the lower right part and the zero matrix to the lower left part.
   for (i = k; i < d; i++)  // If d <= m_order, this does nothing.
      for (j = 0; j < d; j++)
         basis[i][j] = this->m_modulo * (i == j);
   // Fill the rest of the first m_order rows
   for (i = 0; i < dk; i++) {
     for (j = dk; j < d; j++) {
        basis[i][j] = 0;
        for (jj = 1; jj <= m_order; jj++)
           basis[i][j] += m_aCoeff[jj] * basis[i][j - jj] % this->m_modulo;  
     }     
   }
}

//============================================================================

// Builds in basis matrix V^{(p)} directly in dimension 'd' from the entries of the vector 'm_y', 
// as explained in Section 3.1.4 of the guide of LatMRG, puts this matrix in `basis`.  Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0Pol(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   int64_t i, j;
   for (j = 0; j < d; j++) {
      for (i = 0; i < dk; i++)
         basis[i][j] = m_y[j-i+k-1];
   }
   // Fill the rest of the rows
   for (i = dk; i < d; i++) { 
      for (j = 0; j < d; j++)
         basis[i][j] = (i == j) * this->m_modulo;
   }
   
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   this->setDimDual(d);
   this->buildDualBasis0(this->m_dualbasis, d);
   this->setDualNegativeNorm();
}

//============================================================================

// Builds the m-dual basis in a direct way in d dimensions.
// If the basis matrix is in the polynomial form V^{(p)}, then the m-dual basis matrix
// needs to be calculated using mDualUpperTriangular from the class 'BasisConstruction' 
// of 'LatticeTester'. For V^{(0)}, there is a direct way to build the m-dual basis matrix,
// see Sections 3.1.4 and 3.1.5 of the guide.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis0(IntMat &basis, int64_t d) {
   if (use_polynomial_basis) {
      this->buildBasis0Pol(m_copy_primal_basis, d);      
      mDualUpperTriangular(basis, m_copy_primal_basis, this->m_modulo, d);
   }
   else { // Bulids the dual basis according to Eq. (25) in the guide
      this->buildBasis0(m_copy_primal_basis, d);
      assert(d <= this->m_maxDim);
      int64_t k = this->m_order;
      int64_t dk = min(d, k);
      int64_t i, j;
      for (i = 0; i < dk; i++) {
         for (j = 0; j < d; j++)
            basis[i][j] = this->m_modulo * (i == j);  // Avoid "if" statements.
      }
      for (i= dk; i < d; i++) {
         for (j= 0; j <=i; j++)
            basis[i][j] = - m_copy_primal_basis[j][i];
         basis[i][i] = 1;
      }      
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   this->incDimBasis0(this->m_basis, d);
}

//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions, see Section 3.1.6 of the guide.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = m_order;
   int64_t i, j, jj;
   // Add new row to the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d > k) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (jj = 1; jj <= m_order; jj++)
            basis[i][d - 1] += m_aCoeff[jj] * basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
}

//============================================================================

// Increase the dimension of the dual matrix. The algorithm depends on
// whether V^{(0)} or V^{(p)} is used.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDimDual();  // New current dimension. 
   this->setDimDual(d);        
   while (this->m_dim < d) {  // Increase dimension if needed.
      this->m_dim++;
   }
   if (use_polynomial_basis) this->incDimDualBasis0Pol(this->m_dualbasis, d);
   else this->incDimDualBasis0(this->m_dualbasis, d);
}

//============================================================================

// The function increases the dimension of the m-dual basis matrix if we work with
// V^{(0)}, see Section 3.1.5 of the guide.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis0(IntMat &basis, int64_t d) {
   int64_t i;
   incDimBasis0(m_copy_primal_basis, d);
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
   }
   for (i = 0; i < d-1; i++) {
      basis[d-1][i] = -m_copy_primal_basis[i][d-1];  
   } 
   basis[d-1][d-1] = 1;
}

//============================================================================

// If V^{(p)} is used, then m-dual basis matrix needs to be calculated by 
// using mDualUpperTriangular if the order 'k' is greater than 1, see Section 3.1.5
// of the guide.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis0Pol(IntMat &basis, int64_t d) {
   int64_t i;
   incDimBasis0(m_copy_primal_basis, d);
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
   }   
   if (this->m_order == 1) {    
      for (i = 0; i < d-1; i++) {
         basis[d-1][i] = m_copy_primal_basis[i][d-1];  
      }
      basis[d-1][d-1] = 1;
   }
   else { 
      mDualUpperTriangular(m_copy_dual_basis, m_copy_primal_basis, this->m_modulo, d);
      for (i = 0; i < d; i++) {
        basis[d-1][i] = m_copy_dual_basis[d-1][i];  
      }
   }
}

//============================================================================

// We use the columns of basis to construct generating vectors and a pbasis for the projection.
// This function returns the value of `projCase`, which is `true` iff the first m_order coordinates
// are all in the projection, see Section 3.1.7 of the guide. The algorithm depends on which
// basis type (V^{(0) or V^{(p)}) is used.
template<typename Int, typename Real>
bool MRGLattice<Int, Real>::buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis,
      const Coordinates &proj) {
   int64_t d = proj.size();
   int64_t i, j;   
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.

   // Algorithm taylored to the polynomial basis V^{(p)}
   
   if (use_polynomial_basis) {
      int64_t k = this->m_order;
      int64_t dk = min(d, k);
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         // Set column j of all generating vectors, for (j+1)-th coordinate of proj.
         for (i = 0; i < this->m_maxDim; i++) {
            m_genTemp[i][j] = (i == j) * this->m_modulo;
            //if (*it - 1 < (unsigned) d && i <dk)
            //   m_genTemp[i][j] = m_y[*it - 1 - i + k -1];
            if (i < dk) {
               m_genTemp[i][j] = m_y[*it - 1 - i + k -1];
            }
         }
      }
      upperTriangularBasis(pbasis, m_genTemp, this->m_modulo, dimbasis, d);      
   }
   else {
      // Algorithm which does not necessarily use the special form of V^{(p)} but applies
      // for an arbitrary choice of a basis, in particular also for V^{(0)}.
      // Check if we are in case 1.
      // This assumes that the coordinates of each projection are always in increasing order!  ***   
      if (d < (unsigned) m_order) projCase1 = false;
      else {
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++) {
            if (j < m_order) {
               if (*it != unsigned(j + 1)) projCase1 = false;
            } else break;
         }
      }
      if (projCase1) {
         // We first compute the first m_order rows of the projection basis.
         for (i = 0; i < m_order; i++) {
            j = 0;
            for (auto it = proj.begin(); it != proj.end(); it++, j++)
               pbasis[i][j] = basis[i][*it - 1];
         }
         // Then the other rows.
         for (i = m_order; i < d; i++)
            for (j = 0; j < d; j++)
               pbasis[i][j] = this->m_modulo * (i == j);
      } else {
         // In this case we need to use the more general algorithm.
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++) {
            // Set column j of all generating vectors, for (j+1)-th coordinate of proj.
            // for (i = 0; i < dimbasis; i++)
            //   m_genTemp[i][j] = basis[i][*it - 1];
            for (i = 0; i < this->m_maxDim; i++)
               m_genTemp[i][j] = basis[i][*it - 1];
       }
       // std::cout << " Generating vectors: \n" << m_genTemp << "\n";
       upperTriangularBasis(pbasis, m_genTemp, this->m_modulo, dimbasis, d);
      }
   }
   return projCase1;
}

//============================================================================

// The number of generating vectors will be this->m_dim.
// The dimension of the projection will be equal to the cardinality of `proj`
// (see the last statement).
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // int64_t d = proj.size();
   assert(*proj.end() <= uint64_t(this->m_dim));
   // IntMat &pbasis = projLattice.getBasis(); // We want to construct this basis of the projection.
   projLattice.setDim(proj.size());
   this->buildProjection0(this->m_basis, this->m_dim, projLattice.getBasis(), proj);
}

//============================================================================

// The number of generating vectors here will be m_dim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t d = proj.size();     // The dimension of this projection.
   // If dimension is not large enough, we increase it.
   while (*proj.end() > uint64_t(this->m_dim)) {
      this->m_dim++;
   }
   IntMat &pbasis = projLattice.getBasis();  // Reference to basis of projection.
   IntMat &pdualBasis = projLattice.getDualBasis();   // And its m-dual.
   projLattice.setDim(d);
   projLattice.setDimDual(d);

   // We first build a basis for the primal basis of the projection
   this->buildProjection0(this->m_basis, this->m_dim, pbasis, proj);
   // Then we simply calculate its dual.
   mDualUpperTriangular(pdualBasis, pbasis, this->m_modulo, d);
}

//============================================================================

// This function applies phi inverse as described in Section 3.1.2 of the guide, and reverses the coordinates.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol) {
   int i, j, k;
   k = this->m_order;
   IntVec c;
   c.SetLength(k);
   col.SetLength(k);
   for (j = 1; j < k+1; j++) {
      c[j-1] = 0;
      NTL::conv(c[j-1], coeff(rep(pcol), k-j));
      for (i = 1; i < j; i ++) {
         c[j-1] += this->m_aCoeff[i]*c[j - 1 - i];
      }
      c[j-1] = c[j-1] % this->m_modulo;
   }
   for (j = 0; j < k; j++) {
      col[j] = c[k-j-1];
   }
}


//============================================================================
template<typename Int, typename Real>
std::string MRGLattice<Int, Real>::toStringCoef() const {
   return toString(m_aCoeff, 1, this->getDim() + 1);
}

//============================================================================

template class MRGLattice<std::int64_t, double> ;
template class MRGLattice<NTL::ZZ, double> ;
template class MRGLattice<NTL::ZZ, xdouble> ;
template class MRGLattice<NTL::ZZ, quad_float> ;
template class MRGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatMRG

#endif

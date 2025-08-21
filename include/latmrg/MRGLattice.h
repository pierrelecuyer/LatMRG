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
    * This initial primal basis will be upper triangular.
    */
   void buildBasis(int64_t dim);

   /**
    * Builds the m-dual lower triangular basis directly in `dim` dimensions. This `dim` must 
    * not exceed `maxDim`. 
    * In order to build the m-dual basis, the primal basis is needed, even if it has not been
    * built yet. Therefore, a copy of the basis matrix is always created upon building the 
    * m-dual basis. It is stored in the variable 'm_bV0'.
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
    * It requires to increase the dimension of 'm_bV0' as well.
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
   
   virtual void buildDualBasis0(IntMat &basis, int64_t d);      

   virtual void incDimBasis0(IntMat &basis, int64_t d);
   
   virtual void incDimDualBasis0(IntMat &basis, int64_t d);   

   virtual bool buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis, const Coordinates &proj);

   /**
    * Builds the vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used to build the matrix V^{(p)},
    * see Section 3.1.2 of the guide for details.
    */
   void buildy(int64_t dim);  
      
   // Order of this MRG.
   int m_order;

   /**
   * The coefficients \f$a_1, ..., a_k\f$ of the MRG recurrence, a_j stored in `m_aCoeff[j]`.
   */
   IntVec m_aCoeff; 

   /**
    * A vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used in the matrix V^{(p)}.
    */
   IntVec m_y;

   /**
   * For generating the dual basis or increasing its dimension, we need a copy of the
   * the first 'm_order' rows of the primal basis
   */
   IntMat m_bV0;


   /**
    * This auxiliary matrix is used to store the generating vectors of a projections
    * before reducing them into a triangular basis.
    */
   IntMat m_genTemp;
   




   


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
   m_bV0.SetDims(this->m_order, maxDim);   
   m_y.SetLength(maxDim + m_order - 1);
   buildy(maxDim + m_order - 1);
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
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.getNormType()) {
   // this->copy(lat); CW: copy constructor currently not available
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
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

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension isset to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   this->setDim(d);
   this->buildBasis0(this->m_basis, d);
   this->setNegativeNorm();
}

//============================================================================

// Builds the upper-triangular basis V^{(0)} directly in `d` dimensions, as explained in Section 3.1.4 of
// the guide of LatMRG, puts this matrix in `basis`.  Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
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
   /*
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
   */
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
// For V^{(0)}, there is a direct way to build the m-dual basis matrix W^{(0)},
// see Section 3.1.4 of the guide. Moreover, the function stores a copy of the
// first m_order rows of the matrix in the variable m_bV0. This
// is required for being able to increase the dimension
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis0(IntMat &basis, int64_t d) {
   IntMat m_copy_primal_basis; // It might be better to make this a global object! Then we do not need to reinitialize it.
   assert(d <= this->m_maxDim);
   int64_t k = this->m_order;
   int64_t dk = min(k, d);
   int64_t i, j, jj;
   if (k == 1) {
      for (i = 0; i < d; i++)
         basis[0][i] = (i==0) * this->m_modulo;
      for (j = 1; j< d; j++) 
         basis[j][0] = - m_y[j];
      for (i = 1; i < d; i++) {
         for (j=1; j < d; j++) 
            basis[i][j] = (i==j);
      }
   }
   else {
      m_copy_primal_basis.SetDims(this->m_maxDim, this->m_maxDim);
      this->buildBasis0(m_copy_primal_basis, d);      
      mDualUpperTriangular(basis, m_copy_primal_basis, this->m_modulo, d);
   }
   // Moreover, we need to store the first k rows of V^{(0)} in dimension d
   for (i = 0; i < this->m_order; i++)  // If d <= m_order, this does nothing.
       for (j = 0; j < this->m_order; j++)
          m_bV0[j][i] = (i == j);
   for (i = 0; i < dk; i++) {
       for (j = dk; j < d; j++) {
          m_bV0[i][j] = 0;
          for (jj = 1; jj <= m_order; jj++)
             m_bV0[i][j] += m_aCoeff[jj] * m_bV0[i][j - jj] % this->m_modulo;
       }       
    }       
   /* 
   // Bulids the dual basis according to Eq. (25) in the guide
    assert(d <= this->m_maxDim);
    int64_t k = this->m_order;
    int64_t dk = min(d, k);
    int64_t i, j, jj;
    for (i = 0; i < d; i++)  // If d <= m_order, this does nothing.
       for (j = 0; j < d; j++)
          basis[j][i] = (i == j);
    for (i = 0; i < dk; i++) {
       for (j = dk; j < d; j++) {
          basis[j][i] = 0;
          for (jj = 1; jj <= m_order; jj++)
             basis[j][i] += m_aCoeff[jj] * basis[j - jj][i] % this->m_modulo;
          m_bV0[i][j] = basis[j][i];
          basis[j][i] = - basis[j][i];  
       }       
    }      
    for (i = 0; i < dk; i++) {
       for (j = 0; j < dk; j++) {
          basis[j][i] = (i == j) * this->m_modulo;  
          m_bV0[i][i] = 1;
       }
    }   
   */
   
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
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= k) {
         for (jj = 1; jj <= k; jj++)
            basis[i][d - 1] += m_aCoeff[jj] * basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
   // Add me_k as new row to the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   basis[d - 1][d - 1] = this->m_modulo;
}

//============================================================================

// Increase the dimension of the dual matrix.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDimDual();  // New current dimension. 
   this->setDimDual(d);        
   while (this->m_dimdual < d) {  // Increase dimension if needed.
      this->m_dimdual++;
   }
   this->incDimDualBasis0(this->m_dualbasis, d);
}

//============================================================================

// The function increases the dimension of the m-dual basis matrix if we work with
// V^{(0)}, see Section 3.1.5 of the guide.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis0(IntMat &basis, int64_t d) {
   int64_t i;

   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
   }
   // Calculate another column of the copy of the primal basis ...
   for (i = 0; i < this->m_order; i++) {
      m_bV0[i][d - 1] = 0;
      for (int jj = 1; jj <= m_order; jj++)
         m_bV0[i][d - 1] += m_aCoeff[jj] * m_bV0[i][d - 1 - jj] % this->m_modulo;
   }
   // ... and use it to add another row to the dual basis.
   for (i = 0; i < this->m_order; i++) {
      basis[d-1][i] = -m_bV0[i][d-1];  
   } 
   for (i = this->m_order; i < d-1; i++)
      basis[d-1][i] = 0;
   basis[d-1][d-1] = 1;
}


//============================================================================

// We use the columns of basis to construct generating vectors and a pbasis for the projection.
// This function returns the value of `projCase`, which is `true` iff the first m_order coordinates
// are all in the projection, see Section 3.1.7 of the guide. 
template<typename Int, typename Real>
bool MRGLattice<Int, Real>::buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis,
      const Coordinates &proj) {
   
   int64_t d = proj.size();
   int64_t i, j;   
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.

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
         for (auto it = proj.begin(); it != proj.end(); it++, j++) {
            pbasis[i][j] = m_y[*it - 1 - i + this->m_order - 1];
         }
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
         for (i = 0; i < this->m_order; i++) {
               m_genTemp[i][j] = m_y[*it - 1 -i + this->m_order - 1];
         }
    }
    // std::cout << " Generating vectors: \n" << m_genTemp << "\n";
    upperTriangularBasis(pbasis, m_genTemp, this->m_modulo, dimbasis, d);
    // 
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
   projLattice.setDim(proj.size());
   this->buildProjection0(this->m_basis, this->m_dim, projLattice.getBasis(), proj);
}

//============================================================================

// The number of generating vectors here will be m_dim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t d = proj.size();     // The dimension of this projection.
   int64_t i, j;
   // If dimension is not large enough, we increase it.
   while (*proj.end() > uint64_t(this->m_dim)) {
      this->m_dim++;
   }
   projLattice.setDim(d);
   projLattice.setDimDual(d);
   IntMat &pbasis = projLattice.getBasis();  // Reference to basis of projection.
   IntMat &pdualBasis = projLattice.getDualBasis();   // And its m-dual.
   // In the special case where i_1 = k = 1, we get the m-dual
   // basis directly from (30), so there is not need to build the primal basis.
   if (this->m_order==1 && *proj.begin() == 1) {
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         pdualBasis[j][0] = -m_y[*it-1];
      }
      for (j = 0; j < d; j++)
         pdualBasis[i][j] = (i==j) * this->m_modulo;
      for (i = 1; i < d; i++) {
         for (j = 1; j < d; j++)
            pdualBasis[i][j] = (i==j);
      }
   }
   else {
      // Otherwise, we first build a basis for the primal basis of the projection
      this->buildProjection0(this->m_basis, this->m_dim, pbasis, proj);
      // Then we simply calculate its dual.
      mDualUpperTriangular(pdualBasis, pbasis, this->m_modulo, d);
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildy(int64_t dim) {
   int64_t k = m_order;
   int64_t j, jj;
  
   // Builde the vector y according to the algorithm described in Section 3.1.2 of the guide.
   // The vector is stored in the variable 'm_y'.
   for (j = 0; j < dim; j++)
     m_y[j] = 0;
   m_y[k-1] = 1;
   for (j = k; j < dim; j++) {
      for (jj = 1; jj <= m_order; jj++)
         m_y[j] += m_aCoeff[jj] * m_y[j-jj] % this->m_modulo;
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

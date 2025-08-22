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
    * The parameter 'maxdimy' is the maximal length of the vector 'y' which is used
    * to build the basis and the projections. If not set by the user it is equal to maxdim + a.length() - 2
    * so that a basis of dimension 'maxdimy' can be beuilt. As the projections can be built independent of 
    * the basis, the user may want to set a different value of 'maxdimy'.
    * This constructor does not build the basis, so we can build it for a smaller
    * number of dimensions or only for selected projections.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, int64_t maxdimy = 0, NormType norm = L2NORM);

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
   virtual void buildBasis(int64_t dim);

   /**
    * Builds the m-dual lower triangular basis directly in `dim` dimensions. This `dim` must 
    * not exceed `maxDim`. 
    */
   virtual void buildDualBasis(int64_t dim); 

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `maxDim`.
    */
   virtual void incDimBasis();

   /**
    * Increases the current dimension of the m-dual basis by 1.
    * The new increased dimension must not exceed `maxDim`.
    * This function uses the simplified method for MRG lattices given in the lattice tester guide.
    * It requires to increase the dimension of 'm_bV0' as well which is initialized upon the
    * first call of incDimDualBasis. 
    */
   virtual void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. A basis for the projection of this
    * `MRGLattice` over the coordinates in `proj` is returned in `projLattice`.
    * If the first m_order coordinates are in the projection, the algorithm simplifies.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Overrides the same function from `IntLattice`.
    * A lower-triangular basis for the m-dual. If m_order =1 and the first coordinate is included
    * in the projection, things simplify.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /*
    * Returns the dimension up to which the vector 'm_y' has been built.
   */
   int64_t getDimy() {return dim_m_y;};

   /**
    * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
    * where `dim` is the current lattice dimension.
    */
   std::string toStringCoef() const;
   
protected:  

   /**
    * Builds the vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used to build the matrix V^{(p)},
    * see Section 3.1.2 of the guide for details.
    */
   void buildy(int64_t dim);  
   
   /**
    * This function initalizies m_bV0 up to dimension 'dim' 
    */
   void build_m_bV0(int64_t dim);
      
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
    * Current dimension of the vector m_y. It is necessary to store this information 
    * in case we want to build many small dimensional projections but not want to reserve
    * memory for the entire basis matrix or m-dual basis matrix.
    */
   int64_t dim_m_y;

   /**
    * Maximal dimension of m_y
    */
   int64_t m_maxdim_m_y;

   /**
   * For generating the dual basis or increasing its dimension, we need a copy of the
   * the first 'm_order' rows of the primal basis
   */
   IntMat m_bV0;

   /**
    * Stores the current dimension of m_bV0
   */
  int64_t dim_m_bV0;

   /**
    * This auxiliary matrix is used to store the generating vectors of a projections
    * before reducing them into a triangular basis.
    */
   IntMat m_genTemp;

};

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, int64_t maxdimy, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_modulo = m;
   this->m_maxDim = maxDim;   
   setaa(aa);
   this->m_dim = 0;
   m_genTemp.SetDims(maxDim, maxDim); 
   m_bV0.SetDims(this->m_order, maxDim); 
   if (maxdimy == 0)
      maxdimy = maxDim + m_order - 1;
   m_maxdim_m_y = maxdimy;
   m_y.SetLength(maxdimy);
   //buildy(maxDim + m_order - 1);
   buildy(m_order);
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
   // this->copy(lat); // copy constructor currently not available
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_maxdim_m_y = lat.m_maxdim_m_y;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.getNormType()) {
   // this->copy(lat); // copy constructor currently not available
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_maxdim_m_y = lat.m_maxdim_m_y;
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
   this->dim_m_y = 0;
   this->dim_m_bV0 = 0;
}


//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension is set to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   // Make sure that the vector m_y is large enough
   if (getDimy() < d + m_order - 1)
      buildy(d + m_order - 1);
   this->setDim(d);
   assert(d <= this->m_maxDim);
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   int64_t i, j;
   for (j = 0; j < d; j++) {
      for (i = 0; i < dk; i++)
         this->m_basis[i][j] = m_y[j-i+k-1];
   }
   // Fill the rest of the rows
   for (i = dk; i < d; i++) { 
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = (i == j) * this->m_modulo;
   }
   this->setNegativeNorm();
}    

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   this->setDimDual(d);
   assert(d <= this->m_maxDim);
   // Make sure that the vector m_y is large enough
   if (getDimy() < d + m_order - 1) 
      buildy(d + m_order - 1);
   int64_t k = this->m_order;
   int64_t i, j;
   if (k == 1) {
      for (i = 0; i < d; i++)
         this->m_dualbasis[0][i] = (i==0) * this->m_modulo;
      for (j = 1; j< d; j++) 
         this->m_dualbasis[j][0] = - m_y[j] % this->m_modulo;
      for (i = 1; i < d; i++) {
         for (j=1; j < d; j++) 
            this->m_dualbasis[i][j] = (i==j);
      }
   }
   else {
      buildBasis(d);      
      mDualUpperTriangular(this->m_dualbasis, this->m_basis, this->m_modulo, d);
   } 
   this->setDualNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   assert(d <= this->m_maxDim);
   int64_t k = m_order;
   int64_t i, j, jj;
   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = 0;
      if (d - 1 >= k) {
         for (jj = 1; jj <= k; jj++)
            this->m_basis[i][d - 1] += m_aCoeff[jj] * this->m_basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
   // Add me_k as new row to the primal basis.
   for (j = 0; j < d - 1; j++)
      this->m_basis[d - 1][j] = 0;
   this->m_basis[d - 1][d - 1] = this->m_modulo;
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
   int64_t i;
   // Check if m_bV0 already exists and has right dimension
   if (dim_m_bV0 < d - 1)
      build_m_bV0(d-1);
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      this->m_dualbasis[i][d - 1] = 0;
   }
   // Calculate another column of m_bV0 ...
   for (i = 0; i < this->m_order; i++) {
      m_bV0[i][d - 1] = 0;
      for (int jj = 1; jj <= m_order; jj++)
         m_bV0[i][d - 1] += m_aCoeff[jj] * m_bV0[i][d - 1 - jj] % this->m_modulo;
   }
   // ... and use it to add another row to the dual basis.
   for (i = 0; i < this->m_order; i++) {
      this->m_dualbasis[d-1][i] = -m_bV0[i][d-1];  
   } 
   for (i = this->m_order; i < d-1; i++)
      this->m_dualbasis[d-1][i] = 0;
   this->m_dualbasis[d-1][d-1] = 1;
}

//============================================================================

// The number of generating vectors will be this->m_dim.
// The dimension of the projection will be equal to the cardinality of `proj`
// (see the last statement).
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // int64_t d = proj.size();
   // assert(*proj.end() <= uint64_t(this->m_dim));
   projLattice.setDim(proj.size());
   IntMat &pbasis = projLattice.getBasis();  // Reference to basis of projection.
   int64_t d = proj.size();
   int64_t i, j;   
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.
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
   // Make sure that the vector m_y is large enough   
   if (getDimy() < int64_t(*proj.rbegin() + this->m_order - 1)) {
      buildy(*proj.rbegin() + this->m_order - 1);
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
    upperTriangularBasis(pbasis, m_genTemp, this->m_modulo, this->m_dim, d);  
   }   
   
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
   // basis directly from Eq. (30) of the guide, so there is not need to build the primal basis.
   if (this->m_order==1 && *proj.begin() == 1) {
      // Make sure that the vector m_y is large enough
      if (getDimy() < int64_t(*proj.rbegin() + this->m_order - 1)) {
         buildy(*proj.rbegin() + this->m_order - 1);
      }
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
      this->buildProjection(projLattice, proj);
      // Then we simply calculate its dual.
      mDualUpperTriangular(pdualBasis, pbasis, this->m_modulo, d); 
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildy(int64_t dim) {
   int64_t k = m_order;
   int64_t j, jj;  
   // Builds the vector y according to the algorithm described in Section 3.1.2 of the guide.
   // The vector is stored in the variable 'm_y'.
   for (j = 0; j < dim; j++)
     m_y[j] = 0;
   m_y[k-1] = 1;
   for (j = k; j < dim; j++) {
      for (jj = 1; jj <= m_order; jj++)
         m_y[j] += m_aCoeff[jj] * m_y[j-jj] % this->m_modulo;
   } 
   dim_m_y = dim;
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::build_m_bV0(int64_t d) {
   int64_t i, j, jj;
   int64_t dk = min(d, this->m_order);
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
    dim_m_bV0 = d;   
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

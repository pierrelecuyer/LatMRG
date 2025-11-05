#ifndef LATMRG_MWCLATTICE_H
#define LATMRG_MWCLATTICE_H

#include <cassert>
// #include "latticetester/FlexTypes.h"
// #include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
// #include "latticetester/BasisConstruction.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/MWCComponent.h"

namespace LatMRG {

/**
 * This subclass represents a special case of an MRG for which the order `k` is `.
 * Some of the functions can be simplified and/or more efficient in this case.
 * The state of the LCG evolves as \f$$x_n = a x_{n-1} \bmod m\f$$,
 * for a modulus \f$m\f$ and a multiplier \f$\a\f$.
 * The lattices in this class are also special cases of the `Rank1Lattice` objects
 * in Lattice Tester.
 */
template<typename Int, typename Real>
class MWCLattice: public LCGLattice<Int, Real> {

public:

   /**
    * This constructor takes as input the modulus `b`, the vector of multipliers `aa`,
    * the maximal dimension allowed for the full lattice, for a projection, and for a coordinate index,
    * and the norm used to measure the vector lengths.
    * The parameters `maxDim`,  `maxDimProj`, and`maxCoord` are the same as in `LCGLattice`.
    * They give the maximal dimension of a basis (primal or m-dual),
    * the maximal dimension of a projection (for `buildProjection`),
    * and the maximal coordinate index that we may consider for the lattice and for any projection.
    * These values are used by the constructor to allocate space for matrix bases and for the vector
    * \f$\mathbf{y}\f$, which is recomputed each time we set of change the vector `aa`.
    * The constructors do not build a basis.
    */
   MWCLattice(const Int &b, const IntVec &aa, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord,
         NormType norm = L2NORM);

   /**
    * Same as the first constructor, except that `aa` is not given.  It must be set afterwards via `setaa`.
    */
   MWCLattice(const Int &b, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord, NormType norm =
         L2NORM);

   MWCLattice(const Int &b, int64_t maxDim, NormType norm = L2NORM);

   MWCLattice(int64_t maxDim, int64_t maxDimProj, int64_t maxCoord, NormType norm = L2NORM);

   MWCLattice(int64_t maxDim, NormType norm = L2NORM);

   /**
    * Destructor.
    */
   ~MWCLattice();

   /**
    * Sets the vector of multipliers to `aa`.
    * If `maxCoord` is given, `m_maxCoord` is changed to this new value,
    * otherwise its current value is taken.
    */
   virtual void setaa(const IntVec &aa, int64_t maxCoord);

   virtual void setaa(const IntVec &aa);

   /**
    * Sets the value of vector of `b`.
    */
   virtual void setb(const Int &b) {
      m_b = b;
   }

   /**
    * Builds an upper-triangular primal basis in `dim` dimensions, for `dim > 0`
    * and no larger than `m_maxDim` or `m_maxCoord`.
    */
   virtual void buildBasis(int64_t dim) override;

   /**
    * Builds a lower-triangular m-dual basis directly in `dim` dimensions.
    * This `dim` must not exceed `m_maxDim` or `m_maxCoord`.
    */
   virtual void buildDualBasis(int64_t dim) override;

   /**
    * Increases the current dimension of the primal basis by 1 and updates this basis.
    * The new increased dimension must not exceed `m_maxDim` or `m_maxCoord`.
    */
   virtual void incDimBasis() override;

   /**
    * Increases the current dimension of the m-dual basis by 1 and updates this basis.
    * The new increased dimension must not exceed `m_maxDim` or `m_maxCoord`.
    */
   virtual void incDimDualBasis() override;

   /**
    * Builds a basis for the projection of this `MRGLattice` onto the coordinates
    * in `coordSet` and puts it as the `m_basis` of `projLattice`.
    * The construction uses the vector \f$\mathbf{y}\f$, as explained in the guide.
    * The largest coordinate of the projection must not exceed `m_maxCoord` and
    * the number of coordinates must not exceed `m_maxDimProj`.
    * The first coordinate of the projection must be coordinate 1.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` of `projLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

protected:

   Int m_b;  // The MWC modulus.

   // Int m_aa;  // The vector of multipliers, defined in `MRGLattice.h`.

};

//============================================================================
// IMPLEMENTTION

//============================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::MWCLattice(int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      LatMRG::LCGLattice<Int, Real>::LCGLattice(maxDim, maxDimProj, maxCoord, norm) {
}

//============================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::MWCLattice(int64_t maxDim, NormType norm) :
      LatMRG::LCGLattice<Int, Real>::LCGLattice(maxDim, maxDim, maxDim, norm) {
}

//===========================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::MWCLattice(const Int &b, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      MWCLattice<Int, Real>(maxDim, maxDimProj, maxCoord, norm) {
   m_b = b;
}

//===========================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::MWCLattice(const Int &b, int64_t maxDim, NormType norm) :
      MWCLattice<Int, Real>::MWCLattice(b, maxDim, maxDim, maxDim, norm) {
}

//============================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::MWCLattice(const Int &b, const IntVec &aa, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      MWCLattice<Int, Real>(b, maxDim, maxDimProj, maxCoord, norm) {
   setaa(aa);
}

//===========================================================================

template<typename Int, typename Real>
MWCLattice<Int, Real>::~MWCLattice() {
}

//============================================================================
template<typename Int, typename Real>
void MWCLattice<Int, Real>::setaa(const IntVec &aa, int64_t maxCoord) {
   this->m_maxCoord = maxCoord;
   this->m_y.SetLength(maxCoord);
   setaa(aa);
}

template<typename Int, typename Real>
void MWCLattice<Int, Real>::setaa(const IntVec &aa) {
   this->m_aa = aa;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   this->m_order = aa.length()-1;
   this->m_modulo = computeLCGModulusMWC(m_b, aa);
   this->computey(this->m_maxCoord);
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
template<typename Int, typename Real>
void MWCLattice<Int, Real>::buildBasis(int64_t d) {
   assert((d > 0) && (d <= this->m_maxDim) && (d <= this->m_maxCoord));
   this->m_dim = d;
   int64_t i, j;  // Row i and column j.
   for (j = 0; j < d; j++)
      this->m_basis[0][j] = this->m_y[j];
   for (i = 1; i < d; i++)
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = (i == j) * this->m_modulo;
}

//============================================================================

template<typename Int, typename Real>
void MWCLattice<Int, Real>::buildDualBasis(int64_t d) {
   // std::cout << "Start buildDualBasis, d = " << d << ",  maxdim = " << this->m_maxDim << "\n";
   assert((d > 0) && (d <= this->m_maxDim) && (d <= this->m_maxCoord));
   this->setDimDual(d);
   int64_t i, j;
   // Start by putting identity
   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
   // Then make corrections.
   for (j = 0; j < min(d, this->m_order+1); j++)
      this->m_dualbasis[0][j] = this->m_aa[j];
   for (i = 1; i < d; i++)
      this->m_dualbasis[i][i-1] = -m_b;
}

//============================================================================

template<typename Int, typename Real>
void MWCLattice<Int, Real>::incDimBasis() {
   assert((this->m_dim < min(this->m_maxDim, this->m_maxCoord)));
   int64_t d = this->getDim();  // The current (old) dimension.
   this->setDim(d + 1);
   int64_t i, j;
   // Compute the entries of the new column by using the recurrence, if d > 0.
   for (i = 0; i < min(d, 0); i++)
      this->m_basis[i][d] = this->m_a * this->m_basis[i][d - 1] % this->m_modulo;
   // Add `m e_k` as new row to the basis.
   for (j = 0; j < d; j++)
      this->m_basis[d][j] = 0;
   this->m_basis[d][d] = this->m_modulo;
}

//============================================================================

template<typename Int, typename Real>
void MWCLattice<Int, Real>::incDimDualBasis() {
   assert((this->m_dimdual < min(this->m_maxDim, this->m_maxCoord)));
   int64_t d = this->getDimDual();  // The current (old) dimension.
   // std::cout << "Start inDimDualBasis, d = " << d << "\n";
   this->setDimDual(d + 1);
   int64_t i;
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d; i++)
      this->m_dualbasis[i][d] = 0;
   // Then add one more row.
   for (i = 0; i < d-1; i++)
      this->m_dualbasis[d][i] = 0;
   this->m_dualbasis[d][d-1] = -m_b;
   this->m_dualbasis[d][d] = 1;
}

//============================================================================

template<typename Int, typename Real>
void MWCLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   int64_t s = coordSet.size();
   assert(s <= this->m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(this->m_maxCoord));  // Vector y is large enough.
   projLattice.setDim(s);
   IntMat & pbasis = projLattice.getBasis();  // Reference to basis of projection.
   int64_t i, j = 0;
   // Put the first row of the projection basis.
   for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++)
      pbasis[0][j] = this->m_y[*it - 1];
   // Then the other rows.
   for (i = 1; i < s; i++)
      for (j = 0; j < s; j++)
         pbasis[i][j] = this->m_modulo * (i == j);
}

//============================================================================

template<typename Int, typename Real>
void MWCLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   int64_t s = coordSet.size();     // The dimension of this projection.
   assert(s <= this->m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(this->m_maxCoord));  // Vector y is large enough.
   projLattice.setDimDual(s);
   IntMat & pdualBasis = projLattice.getDualBasis();   // Ref to m-dual basis.
   int64_t i, j = 0;
   // First column. The first entry will be overwritten afterwards.
   for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++)
      pdualBasis[j][0] = -this->m_y[*it - 1];
   pdualBasis[0][0] = this->m_modulo;
   // Other columns.
   for (i = 0; i < s; i++)
      for (j = 1; j < s; j++)
         pdualBasis[i][j] = (i == j);
}


} // End namespace LatMRG
#endif

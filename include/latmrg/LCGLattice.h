#ifndef LATMRG_LCGLATTICE_H
#define LATMRG_LCGLATTICE_H

#include <cassert>
// #include "latticetester/FlexTypes.h"
// #include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
// #include "latticetester/BasisConstruction.h"
#include "latmrg/MRGLattice.h"

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
class LCGLattice: public MRGLattice<Int, Real> {

public:

   /**
    * This constructor takes as input the modulus `m`, the multipliers `a`,
    * the maximal dimension allowed for the full lattice, for a projection, and for a coordinate index,
    * and the norm used to measure the vector lengths.
    * The parameters `maxDim`,  `maxDimProj`, and`maxCoord` are the same as in `MRGLattice`.
    * They give the maximal dimension of a basis (primal or m-dual),
    * the maximal dimension of a projection (for `buildProjection`),
    * and the maximal coordinate index that we may consider for the lattice and for any projection.
    * These values are used by the constructor to allocate space for matrix bases and for the vector
    * \f$\mathbf{y}\f$, which is recomputed each time we call `seta`.
    * each time we set of change the vector `aa`.  The constructors do not build a basis.
    */
   LCGLattice(const Int &m, const Int &a, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord,
         NormType norm = L2NORM);

   /**
    * This version takes `maxCoord = maxDimProj = maxDim'.
    */
   LCGLattice(const Int &m, const Int &a, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Same as the first constructor, except that `a` is not given.  It must be set afterwards via `seta`.
    */
   LCGLattice(const Int &m, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord, NormType norm =
         L2NORM);

   LCGLattice(const Int &m, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Copy constructor.
    */
   // LCGLattice(const LCGLattice &Lat);
   /**
    * Assigns `Lat` to this object.
    */
   // LCGLattice& operator=(const LCGLattice &Lat);
   /**
    * Destructor.
    */
   ~LCGLattice();

   /**
    * Sets the multiplier to `a`.
    * This function also recomputes the vector \f$\mathbf{y} = (y_0,\dots,y_{t-1})\f$
    * where \f$t=\f$`maxCoord`\f.
    * If `maxCoord` is given, `m_maxCoord` is changed to this new value,
    * otherwise its current value is taken.
    */
   virtual void seta(const Int &a, int64_t maxCoord);

   virtual void seta(const Int &a);

   /**
    * Compute the entries of \f$(y_0, y_1, ..., y_{t-1})\f$, where \f$\t=\f$`maxCoord` is the largest
    * coordinate index that we want to use for the lattice or for a projection.
    */
   virtual void computey(int64_t maxCoord) override;

   virtual void computey() override;

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
    * in `proj` and puts it as the `m_basis` of `projLattice`.
    * The construction uses the vector \f$\mathbf{y}\f$, as explained in the guide.
    * The largest coordinate of the projection must not exceed `m_maxCoord` and
    * the number of coordinates must not exceed `m_maxDimProj`.
    * The first coordinate of the projection must be coordinate 1.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` of `projLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj)
         override;

protected:

   Int m_a;  // The multiplier.

};

//============================================================================
// IMPLEMENTTION

//============================================================================

template<typename Int, typename Real>
LCGLattice<Int, Real>::LCGLattice(const Int &m, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      LatMRG::MRGLattice<Int, Real>::MRGLattice(m, 1, maxDim, maxDimProj, maxCoord, norm) {
}

//===========================================================================

template<typename Int, typename Real>
LCGLattice<Int, Real>::LCGLattice(const Int &m, int64_t maxDim, NormType norm) :
      LCGLattice<Int, Real>::LCGLattice(m, maxDim, maxDim, maxDim, norm) {
}

//============================================================================

template<typename Int, typename Real>
LCGLattice<Int, Real>::LCGLattice(const Int &m, const Int &a, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      LCGLattice<Int, Real>(m, maxDim, maxDimProj, maxCoord, norm) {
   seta(a);
}

//============================================================================

template<typename Int, typename Real>
LCGLattice<Int, Real>::LCGLattice(const Int &m, const Int &a, int64_t maxDim, NormType norm) :
      LCGLattice<Int, Real>(m, a, maxDim, maxDim, maxDim, norm) {
}

//===========================================================================

template<typename Int, typename Real>
LCGLattice<Int, Real>::~LCGLattice() {
}

//=========================================================================
/*
 template<typename Int, typename Real>
 LCGLattice<Int, Real>::LCGLattice(const LCGLattice &lat) :
 LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(lat.m_modulo, 0, lat.getDim(),
 lat.getNorm()) {
 m_a = lat.m_a;
 m_shift = lat.m_shift;
 }

 //=========================================================================

 template<typename Int, typename Real>
 LCGLattice<Int, Real>& LCGLattice<Int, Real>::operator=(const LCGLattice &lat) {
 if (this == &lat) return *this;
 this->m_dim = lat.m_dim;
 this->copyBasis(lat);
 this->m_order = lat.m_order;
 m_a = lat.m_a;
 m_shift = lat.m_shift;
 return *this;
 }
 */

//============================================================================
template<typename Int, typename Real>
void LCGLattice<Int, Real>::seta(const Int &a, int64_t maxCoord) {
   this->m_maxCoord = maxCoord;
   this->m_y.SetLength(maxCoord);
   seta(a);
}

template<typename Int, typename Real>
void LCGLattice<Int, Real>::seta(const Int &a) {
   m_a = a;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   computey(this->m_maxCoord);
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
template<typename Int, typename Real>
void LCGLattice<Int, Real>::buildBasis(int64_t d) {
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
void LCGLattice<Int, Real>::buildDualBasis(int64_t d) {
   assert((d > 0) && (d <= this->m_maxDim) && (d <= this->m_maxCoord));
   this->setDimDual(d);
   int64_t i, j;
   // Fill first column.
   this->m_dualbasis[0][0] = this->m_modulo;
   for (i = 1; i < d; i++)
      this->m_dualbasis[i][0] = -this->m_y[i] % this->m_modulo;
   // Then the other columns.
   for (i = 0; i < d; i++)
      for (j = 1; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::incDimBasis() {
   assert((this->m_dim < min(this->m_maxDim, this->m_maxCoord)));
   int64_t d = this->getDim();  // The current (old) dimension.
   this->setDim(d + 1);
   int64_t i, j;
   // Compute the entries of the new column by using the recurrence, if d > 0.
   for (i = 0; i < min(d, 0); i++)
      this->m_basis[i][d] = m_a * this->m_basis[i][d - 1] % this->m_modulo;
   // Add `m e_k` as new row to the basis.
   for (j = 0; j < d; j++)
      this->m_basis[d][j] = 0;
   this->m_basis[d][d] = this->m_modulo;
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::incDimDualBasis() {
   assert((this->m_dimdual < min(this->m_maxDim, this->m_maxCoord)));
   int64_t d = this->getDimDual();  // The current (old) dimension.
   this->setDimDual(d + 1);
   int64_t i;
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d; i++)
      this->m_dualbasis[i][d] = 0;
   // Then add the last row of W^{(p)}.
   this->m_dualbasis[d][0] = -this->m_y[d];
   for (i = 1; i < d; i++)
      this->m_dualbasis[d][i] = 0;
   this->m_dualbasis[d][d] = 1;
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t s = proj.size();
   assert(s <= this->m_maxDimProj);
   assert(*proj.end() <= uint64_t(this->m_maxCoord));  // Vector y is large enough.
   projLattice.setDim(s);
   IntMat & pbasis = projLattice.getBasis();  // Reference to basis of projection.
   int64_t i, j = 0;
   // Put the first row of the projection basis.
   for (auto it = proj.begin(); it != proj.end(); it++, j++)
      pbasis[0][j] = this->m_y[*it - 1];
   // Then the other rows.
   for (i = 1; i < s; i++)
      for (j = 0; j < s; j++)
         pbasis[i][j] = this->m_modulo * (i == j);
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t s = proj.size();     // The dimension of this projection.
   assert(s <= this->m_maxDimProj);
   assert(*proj.end() <= uint64_t(this->m_maxCoord));  // Vector y is large enough.
   projLattice.setDimDual(s);
   IntMat & pdualBasis = projLattice.getDualBasis();   // Ref to m-dual basis.
   int64_t i, j = 0;
   // First column. The first entry will be overwrtten after.
   for (auto it = proj.begin(); it != proj.end(); it++, j++)
      pdualBasis[j][0] = -this->m_y[*it - 1];
   pdualBasis[0][0] = this->m_modulo;
   // Other columns.
   for (i = 0; i < s; i++)
      for (j = 1; j < s; j++)
         pdualBasis[i][j] = (i == j);
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::computey(int64_t maxCoord) {
   this->m_maxCoord = maxCoord;
   this->m_y.SetLength(maxCoord);
   computey();
}

//============================================================================

template<typename Int, typename Real>
void LCGLattice<Int, Real>::computey() {
   int64_t j;
   this->m_y[0] = 1;
   for (j = 1; j < this->m_maxCoord; j++)
      this->m_y[j] = m_a * this->m_y[j - 1] % this->m_modulo;
}

//===========================================================================

template class LCGLattice<std::int64_t, double> ;
template class LCGLattice<NTL::ZZ, double> ;
template class LCGLattice<NTL::ZZ, xdouble> ;
template class LCGLattice<NTL::ZZ, quad_float> ;
template class LCGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatMRG
#endif

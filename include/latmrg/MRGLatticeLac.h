#ifndef LATMRG_MRGLATTICELAC_H
#define LATMRG_MRGLATTICELAC_H

#include <string>
#include "latticetester/EnumTypes.h"
#include "latticetester/Types.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"
//#include "latticetester/Types.h"
// #include "latticetester/Lacunary.h"
#include "latticetester/MRGLattice.h"
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
class MRGLatticeLac: public LatticeTester::MRGLattice<Int, Real> {

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
   //MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, IntVec &lac,
   //      NormType norm = L2NORM);
   MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, IntVec & lac, NormType norm = L2NORM);
   


   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac(const MRGLatticeLac &Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac& operator=(const MRGLatticeLac &Lat);

   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac(){ };

   /**
    * Sets the lacunary indices for this lattice to `lac`.
    */
   void setLac(const IntVec &lac);
   
   // CW
   void buildBasis0(IntMat &basis, int64_t d) override;
   void incDimBasis0(IntMat &basis, int64_t d) override;
   bool buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis, const Coordinates &proj) override;

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   Int& getLac(int j);

   /**
    * Takes the polynomial `pcol` and returns in `col` the corresponding column in the
    * matrix of generating vectors.
    */
   void polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol); 


protected:

   // void initStates();


   /**
    * The lacunary indices.
    */
   IntVec m_lac;

   /**
    * The characteristic polynomial `P(z)`.  Maybe no need to store it.
    * We can redefine  `setaa` so it computes and set it.
    * We can also compute it and set it in `buildBasis0`.  Probably safer.  ******
    */
   typename FlexModInt<Int>::PolX m_Pz;   // Maybe not needed.


};

//===========================================================================
// Implementation:


// Constructor.
template<typename Int, typename Real>
MRGLatticeLac<Int, Real>::MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim,
      IntVec &lac, NormType norm) : MRGLattice<Int, Real>(m, aa, maxDim, norm) {
      m_lac = lac;
      ZZ_p::init(m);
}

// This applies phi inverse as described in the guide, and reverses the coordinates.
template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol) {
   int i, j, k;
   k = this->m_order;
   Int temp;
   IntVec c;
   c.SetLength(k);
   col.SetLength(k);
   for (j = 1; j < k+1; j++) {
      NTL::conv(c[j-1], coeff(rep(pcol), k-j));
      for (i = 1; i < j; i ++) {
         NTL::conv(temp, coeff(rep(pcol), j-1-i));
         c[j-1] += this->m_aCoeff[i]*temp;
      }
   }
   for (j = 0; j < k; j++) {
      col[j] = c[k-j-1];
   }
}


//============================================================================

// Replaces the corresponding function in `MRGLattice`.

// Builds a basis directly in `d` dimensions, as explained in Section 4.1.9 of
// the guide of LatMRG.  Must have d <= m_maxdim.
// This should be very similar to `buildProjection0`, except that the indices are usually far apart.

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   int64_t i, j, jj;
   IntVec col;
   
   typename FlexModInt<Int>::PolE polDegOne;
   typename FlexModInt<Int>::PolE polPower;
   
   // Set the characteristic polynomial of the recurrence
   for (int64_t j = 1; j < this->m_aCoeff.length(); j++) {
      SetCoeff(m_Pz, this->m_aCoeff.length() - j - 1, to_ZZ_p(this->m_aCoeff[j]));
   }
   SetCoeff(m_Pz, this->m_aCoeff.length() - 1, 1);
   
   ZZ_pE::init(m_Pz);
   
   //std::cout << m_Pz << "\n";
   
   // Auxilliary variables to calculate powers p^n(z)
   std::string str = "[0 1]";
   std::istringstream in (str);
   in >> polDegOne;   
   //std::cout << polDegOne << "\n";
   // Calculate powers p^\mu-1 for \mu in the set of lacunary indices
   for (j= 0; j < m_lac.length(); j++) {
      power (polPower, polDegOne, m_lac[j] - 1);
      std::cout << m_lac[j] << ": " << polPower << "\n";
      polyToColumn(col, polPower);
      std::cout << col << "\n";
   }
   

   // REDO the following completely.                                             ********
  /*
   for (j = 0; j < k-1; j++)
      this->m_y[j] = 0;
   this->m_y[k-1] = 1;
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
   for (j = k; j < d + k - 1; j++) {
      // Calculate y_j = (a_1 y_{j-1} + ... + a_k y_{j-k}) mod m.
      this->m_y[j] = 0;
      for (jj = 1; jj <= k; jj++)
         this->m_y[j] += this->m_aCoeff[jj] * this->m_y[j - jj];
      this->m_y[j] = this->m_y[j] % this->m_modulo;
      for (i = 0; i < min(k, d - j + k - 1); i++) {  // We want i < k and i+j-k+1 < d.
         basis[i][i + j - k + 1] = this->m_y[j];
      }
   }
   */

}


//=========================================x===================================

// Increases the dimension of given basis from d-1 to d dimensions.
// We compute one new column using the polynomial representation.
//
// REDO this.       *****************
//

template<typename Int, typename Real>
void MRGLatticeLac<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   int64_t i, j, k;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d - 1 >= this->m_order) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= this->m_order) {
         for (k = 1; k <= this->m_order; k++)
            basis[i][d - 1] += this->m_aCoeff[k] * basis[i][d - 1 - k] % this->m_modulo;
      }
   }
}



//============================================================================

// We must be able to do this with basis equal to either m_basis or m_basis0.
// We use the columns of basis to construct generating vectors and a pbasis for the projection.
// This function returns the value of `projCase`, which is `true` iff the first m_order coordinates
// are all in the projection.

// This one may still work with no change. Check carefully and adjust if needed.      ************
//

template<typename Int, typename Real>
bool MRGLatticeLac<Int, Real>::buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis,
      const Coordinates &proj) {
   int64_t d = proj.size();
   int64_t i, j;
   // projCase1 = false;
   // Check if we are in case 1.
   // This assumes that the coordinates of each projection are always in increasing order!  ***
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.
   if (d < (unsigned) this->m_order) projCase1 = false;
   else {
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         if (j < this->m_order) {
            if (*it != unsigned(j + 1)) projCase1 = false;
         } else break;
      }
   }
   if (projCase1) {
      // We first compute the first m_order rows of the projection basis.
      for (i = 0; i < this->m_order; i++) {
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++)
            pbasis[i][j] = basis[i][*it - 1];
      }
      // Then the other rows.
      for (i = this->m_order; i < d; i++)
         for (j = 0; j < d; j++)
            pbasis[i][j] = this->m_modulo * (i == j);
   } else {
      // In this case we need to use the more general algorithm.
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         // Set column j of all generating vectors, for (j+1)-th coordinate of proj.
         for (i = 0; i < dimbasis; i++)
            this->m_genTemp[i][j] = basis[i][*it - 1];
      }
      // std::cout << " Generating vectors: \n" << m_genTemp << "\n";
      upperTriangularBasis(this->m_genTemp, pbasis, this->m_modulo, dimbasis, d);
   }
   return projCase1;
}


// Check if any other function in the parent `MRGLattice`  needs to be changed.      ************

//============================================================================

// The old implementation.  Will be removed.                *********
//
/*
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildLaBasis(int64_t d) {

   if (this->m_order > ORDERMAX)
      LatticeTester::MyExit(1, "MRGLattice::buildLaBasis:   k > ORDERMAX");

   initStates();
   int64_t IMax = m_lac.getSize();

   IntVec b;
   b.SetLength(this->m_order + 1);
   LatticeTester::Invert(m_aCoeff, b, this->m_order);

   // b is the characteristic polynomial of the MRG.
   PolyPE < Int > ::setModulus(this->m_modulo);
   PolyPE < Int > ::setPoly(b);
   PolyPE<Int> pol;
   int64_t ord = 0;

   // Construction d'un systeme generateur modulo m.
   for (int64_t k = 0; k < IMax; k++) {
      // pour chaque indice lacunaire
      NTL::conv(m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod(m_e);
      pol.toVector(m_xi);

      ord = 0;
      for (int64_t i = 1; i <= this->m_order; i++) {
         if (m_ip[i]) {
            ++ord;
            m_t5 = 0;
            for (int64_t j = 1; j <= this->m_order; j++)
               m_t5 += m_sta[i][j] * m_xi[j - 1];
            this->m_wSI[ord][k] = m_t5;
         }
      }
   }

   //  From here we can use BasisConstruction.  *********
*/

   /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
    * les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
    * V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
    * avec m sur la diagonale).
    * */
/*
   LatticeTester::Triangularization<IntMat>(this->m_wSI, this->m_vSI, ord, IMax, this->m_modulo);
   LatticeTester::CalcDual < IntMat > (this->m_vSI, this->m_wSI, IMax, this->m_modulo);

   // Construire la base de dimension 1
   this->m_basis[0][0] = this->m_vSI[0][0];
   this->m_dualbasis[0][0] = this->m_wSI[0][0];
   this->setDim(1);

   this->setNegativeNorm();
   this->setDualNegativeNorm();

   //  This approach could be slow!  *****
   for (int64_t i = 2; i <= d; i++)
      incDimLaBasis(IMax);

   // for debugging
   // trace("ESPION_1", 1);
}
*/

//============================================================================
/*
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimLaBasis(int64_t IMax) {

   LatticeTester::IntLatticeExt<Int, Real>::incDim();
   const int64_t dim = this->getDim(); // new dimension (dim++)
*/
   /*
    if (dim >= IMax) {
    MyExit (0,
    "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
    }
    */
/*
   IntVec tempLineBasis(dim);
   IntVec tempColBasis(dim);

   for (int64_t i = 0; i < dim - 1; i++) {

      // tempLineBasis <- m_basis[i]
      for (int64_t k = 0; k < dim - 1; k++)
         tempLineBasis[k] = this->m_basis[i][k];

      for (int64_t i1 = 0; i1 < dim - 1; i1++) {

         Int tempScalDual;
         LatticeTester::ProdScal<Int>(tempLineBasis, this->m_wSI[i1], dim, tempScalDual);
         LatticeTester::Quotient(tempScalDual, this->m_modulo, tempScalDual);
         this->m_t1 = tempScalDual * this->m_vSI[i1][dim - 1];
         tempColBasis[i] += this->m_t1;
      }
      LatticeTester::Modulo(tempColBasis[i], this->m_modulo, tempColBasis[i]);
      this->m_basis[i][dim - 1] = tempColBasis[i];
   }

   for (int64_t j = 0; j < dim - 1; j++)
      this->m_basis[dim - 1][j] = 0;
   this->m_basis[dim - 1][dim - 1] = this->m_vSI[dim - 1][dim - 1];

   for (int64_t i = 0; i < dim - 1; i++)
      this->m_dualbasis[i][dim - 1] = 0;

   for (int64_t j = 0; j < dim - 1; j++) {

      Int tempScalDualBis;

      for (int64_t i = 0; i < dim - 1; i++) {
         this->m_t1 = this->m_dualbasis[i][j];
         this->m_t1 *= tempColBasis[i];
         tempScalDualBis += this->m_t1;
      }
      if (tempScalDualBis != 0) tempScalDualBis = -tempScalDualBis;

      LatticeTester::Quotient(tempScalDualBis, this->m_vSI[dim - 1][dim - 1], tempScalDualBis);
      this->m_dualbasis[dim - 1][j] = tempScalDualBis;
   }

   LatticeTester::Quotient(this->m_modulo, this->m_vSI[dim - 1][dim - 1], this->m_t1);
   this->m_dualbasis[dim - 1][dim - 1] = this->m_t1;

   this->setNegativeNorm();
   this->setDualNegativeNorm();

}
*/
}
#endif


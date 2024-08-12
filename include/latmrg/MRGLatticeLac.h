#ifndef LATMRG_MRGLATTICELAC_H
#define LATMRG_MRGLATTICELAC_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/Types.h"
// #include "latticetester/Const.h"
#include "latticetester/Lacunary.h"
#include "latticetester/MRGLattice.h"
// #include "latmrg/MRGComponent.h"
#include <string>

namespace LatMRG {

/**
 * This subclass of `MRGLattice` constructs and handles lattice bases built from MRGs as in `MRGLattice`,
 * but with arbitrary lacunary indices that are regularly spaced by packets of the same size.
 *
 * Perhaps the new functions offered here could be integrated into MRGLattice,
 * to reduce the number of classes.  ???            ************
 *
 */
class MRGLatticeLac: public LatticeTester::MRGLattice {
public:

   /**
    * Constructor with modulus of congruence \f$m\f$, order of the recurrence
    * \f$k\f$, multipliers \f$A\f$, maximal dimension `maxDim`, and lattice type
    * `lattype`. Vector and matrix indices vary from 1 to `maxDim`. The length of
    * the basis vectors is computed with `norm`. The bases are built using the
    * *lacunary indices* `lac`.
    * The basis is built for the lacunary
    * indices `lac`.  `aa` has to be a vector of k+1 components
    * with `a[i]`=\f$a_i\f$ for compatibility with other classes.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, int64_t k,
         IntVec &lac, bool withPrimal = false, bool withDual = false,
         LatticeType lattype, LatticeTester::NormType norm =
         LatticeTester::L2NORM);

   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac(const MRGLatticeLac &Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac& operator=(const MRGLatticeLac &Lat);

   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac();

   /**
    * Builds the basis of the MRG recurrence in \f$d\f$ dimensions, not using
    * the lacunary indices.
    */
   void buildBasis(int d);

   /**
    * Builds the basis of the MRG recurrence in case of lacunary indices.
    */
   void buildLaBasis(int64_t d);

   /**
    * Increments the basis by 1 in case of lacunary indices.
    * Uses the method described in the article: P. L'Ecuyer and R. Couture,
    * "An Implementation of the Lattice and Spectral Tests for Multiple
    * Recursive Linear Random Number Generators", INFORMS Journal on
    * Computing, 9, 2 (1997), page 206--217. Section 3, "Lacunary indices".
    */
   void incDimLaBasis(int64_t);

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   Int& getLac(int j);

   /**
    * Sets the lacunary indices for this lattice to `lat`.
    */
   void setLac(const Lacunary &lat);

protected:

   void initStates();

   /**
    * Increases the dimension of the basis by 1.
    */
   void incDimBasis(int);

   /**
    * The lacunary indices.
    */
   Lacunary m_lac;

   //===========================================================================

   /** Max order for lacunary case in this class; otherwise, it takes too much memory.
    For order > ORDERMAX, use subclass MRGLatticeLac instead.
    This means that we can have short lacunary indices, supported here,
    and also long lacunary indices (e.g., for multiple streams), supported in MRGLatticeLac.  ******* ???
    */
#define ORDERMAX 100

};

//===========================================================================
// Implementation:

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildLaBasis(int64_t d) {

   // NOT USED, see: MRGLatticeLac::buildBasis

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

   /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
    * les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
    * V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
    * avec m sur la diagonale).
    * */
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

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimLaBasis(int64_t IMax) {

   LatticeTester::IntLatticeExt<Int, Real>::incDim();
   const int64_t dim = this->getDim(); // new dimension (dim++)

   /*
    if (dim >= IMax) {
    MyExit (0,
    "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
    }
    */

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

}
#endif


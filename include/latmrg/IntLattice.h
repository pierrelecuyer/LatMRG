#ifndef LATMRG__INTLATTICE_H
#define LATMRG__INTLATTICE_H
#include "latticetester/IntLatticeBasis.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

/**
 * \copydoc LatticeTester::IntLattice
 * Beside, it contains the order k of the genrator.
 * In fact, the basis has a particular form which need
 * the order k.
 */
class IntLattice : public LatticeTester::IntLatticeBasis {
public:

   /**
    * \copydoc LatticeTester::IntLattice::IntLattice(const MScal&, int, int, NormType)
    */
   IntLattice (MScal modulo, int k, int maxDim, LatticeTester::NormType norm = LatticeTester::L2NORM);


   /**
    * \copydoc LatticeTester::IntLattice::IntLattice(const IntLattice&)
    * Erwan
    */
   IntLattice (const IntLattice & Lat);

   /**
    * Init the matrix
    */
   void init ();

   /**
    * Return the order of the matrix
    */
   int getOrder() const { return m_order; }

   /**
    * Increment the dimension of the element of the lattice
    */
   void incrementDimension ();

   /**
    * Computes the logarithm of the normalization factor
    * (<tt>m_lgVolDual2</tt>) in all dimensions \f$\le\f$ `MaxDim` for
    * the lattice. `lgm2` is the logarithm in base 2 of \f$m^2\f$.
    */
   void calcLgVolDual2 (double lgm2);


   /**
    * Destructor.
    */
   virtual ~IntLattice ();

#ifdef WITH_NTL
   /**
    * The components of the lattice when it is built out of more than one
    * component. When there is only one component, it is unused as the
    * parameters are the same as above.
    */
   std::vector<MRGComponent *> comp;
#endif


protected:

   /**
    * \copydoc LatticeTester::IntLattice::kill()
    */
   virtual void kill ();

   /**
    * The order of the basis.
    */
   int m_order;

   /**
    * The maximum Dimension for the test
    */
   //int m_maxDim;

   /*
    * Represente sur dual along the diagonal?? ERWAN
    */
   double *m_lgVolDual2;

   /**
    * The logarithm \f$\log_2 (m^2)\f$.
    */
   double m_lgm2;

   /**
    * Use to save the dual basis and the basis in some works.
    */
   BMat m_wSI, m_vSI;

   /**
    * Working Variables use in MRGLattice.h
    */
   MScal m_t1, m_t2, m_t3;

};

}
#endif

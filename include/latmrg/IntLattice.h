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
    */
   IntLattice (const IntLattice & Lat);


   /**
    * Return the order of the matrix
    */
   int GetOrder() const { return m_order; }

   /**
    * Increment the dimension of the element of the lattice
    */
   void IncrementDimension ();

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

};

}
#endif

#ifndef MRGLATTICELAC_H
#define MRGLATTICELAC_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latmrg/Lacunary.h"
#include "latmrg/Const.h"
#include "latmrg/MRGLattice.h"


namespace LatMRG {

/**
 * This class implements lattice bases built from multiple recursive linear
 * congruential generators (see class <tt>MRGLattice</tt>) using *lacunary
 * indices*.
 *
 */
class MRGLatticeLac:   public MRGLattice {
public:

   /**
    * Constructor with modulus of congruence \f$m\f$, order of the recurrence
    * \f$k\f$, multipliers \f$A\f$, maximal dimension `maxDim`, and lattice type
    * `latt`. Vector and matrix indices vary from 1 to `maxDim`. The length of
    * the basis vectors is computed with `norm`. The bases are built using the
    * *lacunary indices* `lac`.
    */
   MRGLatticeLac (const MScal & m, const MVect & A, int maxDim, int k,
                     BVect & lac, LatticeType latt,
                     LatticeTester::NormType norm = LatticeTester::L2NORM);

   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac (const MRGLatticeLac & Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   MRGLatticeLac & operator= (const MRGLatticeLac & Lat);

   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac();

   /**
    * Builds the basis of the MRG recurrence in dimension \f$d\f$ using
    * the lacunary indices.
    */
   void buildBasis (int d);

   /**
    * Increases the dimension of the basis by 1.
    */
   void incDim() { incDimBasis (getDim()); }

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   BScal & getLac (int j);

   /**
    * Sets the lacunary indices for this lattice to `lat`.
    */
   void setLac (const Lacunary & lat);

protected:

   void initStates();

   /**
    * Increases the dimension of the basis by 1.
    */
   void incDimBasis (int);

   /**
    * The lacunary indices.
    */
    Lacunary m_lac;
};

}
#endif

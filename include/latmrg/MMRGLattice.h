#ifndef MMRGLATTICE_H
#define MMRGLATTICE_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latmrg/Lacunary.h"
#include "latmrg/Const.h"
#include "latmrg/IntLattice.h"
#include <string>


namespace LatMRG {

/**
 * This class implements lattice basis built from M-MRG (matrix multiple recursive 
 * linear congruential generators). One must first call the constructor with a
 * given congruence modulus \f$m\f$, a given generator matrix for the 
 * recurrence, and a maximal dimension for the basis. One must then build the
 * lattice basis associated to the generator matrix for a given dimension.
 * Each MMRG is defined by a generator matrix \f$A\f$. This MMRG satisfies 
 * the recurrence
 * \f[
 *   X_n = A X_{n-1} \mod m.
 * \f]
 */
class MMRGLattice: public LatMRG::IntLattice {
public:

   /**
    * Constructor with modulus of congruence \f$m\f$, generator matrix \f$A\f$, dimension 
    * of generator matrix \f$r\f$, maximal dimension `MaxDim`, and lattice type
    * `Latt`. Vectors and (square) matrices of the basis have maximal dimension
    * `maxDim`, and the indices of vectors and matrices vary from dimension 0 to
    * `maxDim`-1. The norm to be used for the basis vectors is `norm`.
    */
   MMRGLattice (const MScal & m, const MMat & A, int maxDim, int r,
                 LatticeType lat, LatticeTester::NormType norm = LatticeTester::L2NORM);

   /**
    * As in the constructor above but the basis is built for the lacunary
    * indices `lac`.
    */
  //PW_TODO à faire plus tard
  MMRGLattice (const MScal & m, const MMat & A, int maxDim, int r, BVect & lac,
               LatticeType lat, LatticeTester::NormType norm = LatticeTester::L2NORM);

   /**
    * Copy constructor. The maximal dimension of the created basis is set
    * equal to <tt>Lat</tt>’s current dimension.
    */
   MMRGLattice (const MMRGLattice & Lat);

   /**
    * Destructor.
    */
   ~MMRGLattice();

   /**
    * Cleans and releases memory used by this object.
    */
   void kill();

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set equal to <tt>Lat</tt>’s current dimension.
    */
   MMRGLattice & operator= (const MMRGLattice & Lat);

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   BScal & getLac (int j);

   /**
    * Sets the lacunary indices for this lattice to `lat`.
    */
   virtual void setLac (const Lacunary & lat);

   /**
    * Returns the generator matrix \f$A\f$ as a string.
    */
   std::string toStringGeneratorMatrix() const;

   /**
    * Builds the basis in dimension \f$d\f$.
    */
   virtual void buildBasis (int d);

   /**
    * Increments the dimension of the basis by 1 by calling either
    * `incDimBasis` or `incDimLaBasis`.
    */
   virtual void incrementDim();

   /**
    * Returns `true` for the case of lacunary indices, returns `false` for
    * non-lacunary indices.
    */
   bool isLacunary() const { return m_lacunaryFlag; }

   /**
    * Returns a non-mutable copy of the generator matrix of the MMRG
    */
   const MMat & getGeneratorMatrix() const { return m_A; }

protected:

   /**
    * Initializes some of the local variables.
    */
   void init();

   /**
    * Builds the basis of the MMRG recurrence in case of non-lacunary
    * indices.
    */
   void buildNonLacunaryBasis (int d);

   /**
    * Builds the basis of the MMRG recurrence in case of lacunary indices.
    */
   void buildLacunaryBasis (int d);

   /**
    * Increments the basis by 1 in case of non-lacunary indices.
    */ 
   //PW_TODO c'était virtual avant : normal ?
   void incrementDimBasis ();

   /**
    * Increments the basis by 1 in case of lacunary indices.
    */
   void incrementDimLacunaryBasis (int);

   /**
    * The generator matrix of the recurrence.
    */
   MMat m_A;

   /**
    * Indicates which lattice or sublattice is analyzed.
    */
   LatticeType m_latType;

   /**
    * Is `true` in the case of lacunary indices, `false` otherwise.
    */
   bool m_lacunaryFlag;

   /**
    * Contains the lacunary indices when `LacunaryFlag` is `true`,
    * otherwise is undefined.
    */
   Lacunary m_lac;




   /**
    * Work variables.
    *
    * @{
    */
    MScal m_t4, m_t5, m_t6, m_t7, m_t8, m_e;
    MVect m_xi;
    /**
     * @}
     */

   /**
    * \f$\clubsuit\f$ Seems to be use as working variables. To be completed. Erwan
    */
   BMat m_sta;

   /**
    * When the flag <tt>m_ip[i]</tt> is `true`, the \f$i\f$-th diagonal
    * element of matrix <tt>m_sta</tt> is non-zero (modulo \f$m\f$) and
    * divides \f$m\f$. Otherwise (when <tt>m_ip[i]</tt> is
    * <tt>false</tt>), the \f$i\f$-th line of matrix <tt>m_sta</tt> is
    * identically 0.
    */
   bool *m_ip;
};

}
#endif

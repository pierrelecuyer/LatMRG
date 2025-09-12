#ifndef	LATMRG_MLCGCOMPONENT_H
#define	LATMRG_MLCGCOMPONENT_H

#include <NTL/mat_poly_ZZ.h>
#include <string>
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"

namespace LatMRG {

using namespace LatMRG;

/**
 * This class offers tools to verify if an MLCG recurrence or the form
 * \f[
 *   \mathbb{x}_n = \mathbb{A} \mathbb{x}_{n-1} \mod m
 * \f]
 * has full period or not. This depends on \f$m\f$ and on the characteristic
 * polynomial of the \f$k\times k\f$ transition matrix \f$\mathbb{A}\f$,
 * and it requires the factorizations of \f$m-1\f$ and of \f$r = (m^k-1)/(m-1)\f$,
 * which can be done once for all if we want to examine several matrices \f$\mathbb{A}\f$.
 */
template<typename Int> class MLCGComponent {

public:

   /**
    * Constructor for modulus \f$m\f$ and \f$k\times k\f$ transition matrices.
    * The matrix \f$\mathbb{A}\f$ must be set separately by `setA`.
    * The Arguments `decom1` and `decor` specify the type of the prime factor
    * decomposition of \f$m-1\f$ and \f$r=(m^k-1)/(m-1)\f$, respectively.
    * The type `DecompType` is defined in `EnumTypes` and offers the
    * following choices:  `DECOMP`: the integer will be factored,
    * `DECOMP_WRITE`: it will be factored and the prime factors written in a file,
    * `DECOMP_READ`: the prime factors are read from a file,
    * `DECOMP_PRIME`: the integer is assumed to be prime.
    * The names `filem1` and `filer` (character strings) specify the file where the factors of
    * \f$m-1\f$ and \f$r\f$, are read or written when the type
    * <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt> is used.
    * The given files must be accessible by the program.
    * The format for the factorization in these files is described in class `IntFactorization`.
    */
   MLCGComponent(const Int &m, int64_t k, DecompType decompm1, const char *filem1, DecompType decompr,
         const char *filer);

   /**
    * Same as the previous constructor, but with \f$m = b^e + c\f$.
    */
   MLCGComponent(Int b, int64_t e, Int c, int64_t k, DecompType decompm1, const char *filem1,
         DecompType decompr, const char *filer);

   /**
    * Destructor.
    */
   ~MLCGComponent();

   /**
    * Sets the matrix \f$\mathbb{A}\f$ to `A`, which must be a \f$k\times k\f$ matrix
    * that contains the elements of \f$\mathbb{A}\f$,
    * If we only want to call `maxPeriod(A)`, we can just pass it \f$\mathbb{A}\f$ and
    * there is no need to call `setA` before.
    */
   void setA(const IntMat &A);

   /**
    * Returns the matrix \f$\mathbb{A}\f$ of the recurrence, in same format as in `setA` above.
    */
   IntMat getA() const {
      return m_A;
   }

   /**
    * Returns `true` if `A` give an MLCG with maximal
    * period; returns `false` otherwise.
    * The internal matrix \f$\mathbb{A}\f$ is also set to this `A`.
    */
   bool maxPeriod(const IntMat &A);

   /**
    * Assumes that the maximal-period condition (1) holds and checks only for
    * conditions 2 and 3. Returns `true` iff these two conditions hold for the
    * matrix \f$\mathbb{A}\f$.
    * The internal matrix \f$\mathbb{A}\f$ is also set to this `A`.
    */
   bool maxPeriod23(const IntMat &A);

   /**
    * Returns the value of the modulus \f$m\f$ of the recurrence.
    */
   const Int getModulus() const {
      return m_m;
   }

   /**
    * Returns the value of \f$k\f$, the dimension of the matrix.
    */
   int64_t getk() const {
      return m_k;
   }

   /**
    * Returns the value of `b` for which `b^e+r = m`, when available.
    */
   Int& getb() {
      return m_b;
   }

   /**
    * Returns the value of `e` for which `b^e+r = m`, when available.
    */
   int64_t& gete() {
      return m_e;
   }

   /**
    * Returns the value of `c` for which `b^e+c = m`, when available.
    */
   Int& getc() {
      return m_c;
   }

   /**
    * Returns a const reference to the initial state `orbitSeed`, when available.
    * */
//   IntVec& getOrbitSeed() {
//      return orbitSeed;
//   }

   /**
    * Returns a descriptor of this object as a string.
    */
   // std::string toString();

private:

   /**
    * This is called by the constructor, with the same arguments.
    */
   void init(const Int &m, int64_t k, DecompType decom1, const char *filem1, DecompType decor,
         const char *filer);

   /**
    * The type of generator this stores. Should be MRG, MMRG or MWC. The
    * default value is LCG.
    */
   // GenType m_type = LCG;   // This should not be here, but should depend on the subclass, it seems.

   /**
    * The prime factor decomposition of \f$m-1\f$.
    */
   IntFactorization<Int> ifm1;

   /**
    * The prime factor decomposition of \f$m\f$.
    * Used to compute the full period of a LCG with a carry.   // Not here!    *******
    */
   // IntFactorization<Int> factor;
   /**
    * The prime factor decomposition of \f$m-1\f$.
    * This is used when computing the full period of a matrix generator.
    * This is because we use NTL to compute the characteristic polynomial
    * of the matrix since it is already implemented.
    */
   // IntFactorization<NTL::ZZ> ifm2;

   /**
    * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
    * \f$k\f$ is the order of the recurrence.
    */
   IntFactorization<Int> ifr;

   /**
    * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
    * \f$k\f$ is the order of the recurrence.
    * This is used when computing the full period of a matrix generator.
    * This is because we use NTL to compute the characteristic polynomial
    * of the matrix since it is already implemented.
    */
   // IntFactorization<NTL::ZZ> ifr2;
   /**
    * The modulus \f$m\f$ of the recurrence, as a  `Modulus` object.
    */
   // Modulus<Int> modulus;    //  Want this ????
   /**
    * The modulus \f$m\f$ of this MRG.
    */
   Int m_m;

   /**
    * The order \f$k\f$ of the recurrence.
    */
   int64_t m_k = 1;

   /**
    * The transition matrix \f$\mathbf{A}\f$.
    */
   IntMat m_A;

   /**
    * Basis `b` with `b^e+c = m`.
    * */
   Int m_b;

   /**
    * Exponent `e` with `b^e+c = m`.
    * */
   int64_t m_e;

   /**
    * Rest `c` with `b^e+c = m`.
    * */
   Int m_c;

   /**
    * The period length \f$\rho\f$ for this MRG. For now, the
    * recurrence is assumed to correspond to a primitive polynomial and
    * `rho` is calculated as
    * \f[
    *   \rho= m^k - 1
    * \f]
    * This value is calculated by `MRGLatticeFactory` and stored here for
    * simplicity.                     Used?   *****************
    */
   Int rho;

   /**
    * Value needed for the calculation of the multipliers of a combined
    * MRG. It is defined by
    * \f[
    *   n_j = (m/m_j)^{-1} \mbox{ mod } m_j \qquad\mbox{ for } j = 1,…,J,
    * \f]
    * where \f$n_j = \f$ `nj`, \f$m_j\f$ is this object’s modulus `m`,
    * \f$m\f$ is the calculated modulus for the combined MRG (see class
    * <tt>MRGLatticeFactory</tt>), and \f$(m/m_j)^{-1} \mbox{ mod } m_j\f$
    * is the inverse of \f$m/m_j\f$ modulo \f$m_j\f$. This value is
    * calculated by `MRGLatticeFactory` and stored here for simplicity.
    */
   Int nj;     // Not sure if we need this here.   ***************

};
// End class declaration

//===========================================================================
// IMPLEMENTTION

//=============================================================================
// Main constructor.
template<typename Int>
MLCGComponent<Int>::MLCGComponent(const Int &m, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   m_b = m;
   m_e = 1;
   m_c = Int(0);
   init(m, k, decompm1, filem1, decompr, filer);
}

//===========================================================================

// Constructor with alternative format for m.
template<typename Int>
MLCGComponent<Int>::MLCGComponent(Int b, int64_t e, Int c, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   Int m = NTL::power(b, e) + c;
   m_b = b;
   m_e = e;
   m_c = c;
   init(m, k, decompm1, filem1, decompr, filer);
}

//============================================================================
template<typename Int>
void MLCGComponent<Int>::init(const Int &m, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   m_m = m;
   setModulusIntP<Int>(m_m);
   m_k = k;
   m_A.SetDims(k, k);
   // orbitSeed.SetLength(m_k);   // Used ???
   Int mm1;  // m-1
   mm1 = m - 1;
   ifm1 = IntFactorization<Int>(mm1);
   ifm1.decompToFactorsInv (decompm1, filem1);
   Int r;
   r = (NTL::power(m, m_k) - 1) / (m - 1);
   ifr = IntFactorization<Int>(r);
   ifr.decompToFactorsInv (decompr, filer);
}

//===========================================================================
template<typename Int>
MLCGComponent<Int>::~MLCGComponent() {
   //a.kill();
   //orbitSeed.kill();
}

//===========================================================================

template<typename Int>
void MLCGComponent<Int>::setA(const IntMat &A) {
   assert(A.NumRows() == m_k);
   m_A = A;
}

//===========================================================================

template<typename Int>
bool MLCGComponent<Int>::maxPeriod(const IntMat &A) {
   setA(A);
   return isPrimitive(A, m_m, ifm1, ifr);
}

//===========================================================================
// Must be rewritten!
template<typename Int>
bool MLCGComponent<Int>::maxPeriod23(const IntMat &A) {
   setA(A);
   return isPrimitive23(A, m_m, ifr);
}

template class MLCGComponent<std::int64_t> ;
template class MLCGComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif

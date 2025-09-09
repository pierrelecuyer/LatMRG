#ifndef	LATMRG_LCGCOMPONENT_H
#define	LATMRG_LCGCOMPONENT_H

#include <string>
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"

namespace LatMRG {

/**
 * This class offers tools to test the period of an LCG recurrence modulo \f$m\f$, of the form
 * \f[
 *   x_n = (a x_{n-1} + c_0) \bmod m.
 * \f]
 * The object is constructed for a given value of \f$m\f$ and a boolean that tells if we assume a nonzero increment
 * \f$c_0\f$ relatively prime with \f$m\f$, or no increment (\f$c_0 = 0\f$).
 * The constructor takes the value of \f$m\f$ and computes or read the required factorizations.
 * The function `maxPeriod` tests the period for a given value of \f$a\f$.
 * This can be done for several values of \f$a\f$ with the same `LCGComponent` object.
 * No value of \f$a\f$ is saved in the object.
 */
template<typename Int>
class LCGComponent {

//  using namespace LatMRG;

public:

   /**
    * Construct an object for modulus \f$m\f$, with the type of prime factor decomposition
    * of \f$m-1\f$ or \f$m\f$ specified by `decomp`.
    * When `increment` is `false`, we consider an LCG with no increment and we need the
    * decomposition of \f$m-1\f$, otherwise we consider an LCG with an increment \f$c_0\f$
    * relatively prime with \f$m\f$ and we need the prime decomposition of \f$m\f$.
    * The type `DecompType` is defined in `EnumTypes` and offers the following choices:
    * `DECOMP`: the integer will be factored,
    * `DECOMP_WRITE`: it will be factored and the prime factors written in a file,
    * `DECOMP_READ`: the prime factors are read from a file,
    * `DECOMP_PRIME`: the integer is assumed to be prime.
    * The name `filename` specifies the file where the factors of
    * \f$m-1\f$ or \f$m\f$ are read or written when the type
    * <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt> is used.
    * The given files must be accessible by the program.
    * The format for the factorization in these files is described in class `IntFactorization`.
    */
   LCGComponent(const Int &m, DecompType decomp, const char *filename, bool increment = false);

   /**
    * Same as the previous constructor, but with `m=b^e+c`.
    */
   LCGComponent(const Int &b, int e, const Int &c, DecompType decomp, const char *filename, bool increment = false);

   /**
    * Destructor.
    */
   ~LCGComponent();

   /**
    * Returns `true` iff the LCG has full period. If there is no increment \f$c_0\f$,
    * the function checks if f$a\f$ is a primitive element modulo \f$m\f$.
    * If there is an increment, this increment is assumed to be
    * relatively prime with \f$m\f$ and this function checks the two conditions:
    * (1) Every prime divisor \f$q\f$ of \f$m\f$ must divide \f$a-1\f$;
    * (2) If 4 divides \f$m\f$, then it must divide \f$a-1\f$.
    */
   bool maxPeriod(const Int &a);

   /**
    * Returns the value of the modulus \f$m\f$.
    */
   const Int getModulus() const {
      return m_m;
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
   int& gete() {
      return m_e;
   }

   /**
    * Returns the value of `c` for which `b^e+c = m`, when available.
    */
   Int& getc() {
      return m_c;
   }

   /**
    * Returns the factorization of `m` or `m-1`.
    */
   IntFactorization<Int>& getFactorization() {
      return m_fact;
   }

   /**
    * Returns a const reference to the initial state `orbitSeed`, when available.
    */
   //Int& getOrbitSeed() {
   //   return orbitSeed;
   //}

   /**
    * Returns this object as a string.
    */
   std::string toString();

   /**
    * The equivalent modulo of the MWC generator, if relevant. This is because
    * we check MWC period with the modulo of the equivalent LCG.                   ?????
    */
   //  Int m_MWCb;

private:

   /**
    * This is called by the constructor, with the same arguments.
    */
   void init(const Int &m, DecompType decomp, const char *filename, bool increment = false);

   /**
    * Used for the prime factor decomposition of \f$m-1\f$ (when there is no increment)
    * or \f$m\f$ (when the LCG has an increment).
    */
   IntFactorization<Int> m_fact;

   /**
    * The modulus \f$m\f$ of the recurrence.
    */
   Int m_m;

   /**
    * Basis `b` when `m = b^e + c`.
    */
   Int m_b;

   /**
    * Exponent `e` when `m = b^e + c`.
    */
   int m_e;

   /**
    * Rest `c` when `m = b^e + c`.
    */
   Int m_c;

   /**
    * Indicates if the LCG has an additive increment or not.
    */
   bool m_increment = false;

   /**
    * The length of the period \f$\rho\f$ for this MRG. For now, the
    * recurrence is assumed to correspond to a primitive polynomial and
    * `rho` is calculated as
    * \f[
    *    \rho = m - 1
    * \f]
    * This value is calculated by `MRGLatticeFactory` and stored here for
    * simplicity.
    */
   // Int rho;
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
   // Int nj;

   /**
    * Contains the starting state of the component for the case when the
    * lattice type is `ORBIT`. It is made of \f$k\f$ numbers.
    */
   // Int orbitSeed;   // Not used for now.  Should be in LCGLattice.

};   // End of class declaration



//===========================================================================
// IMPLEMENTATION

//===========================================================================

template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &m, DecompType decomp, const char *filename, bool increment) {
   m_b = m;
   m_e = 1;
   m_c = Int(0);
   init(m, decomp, filename, increment);
}

//===========================================================================

template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &b, int e, const Int &c, DecompType decomp, const char *filename, bool increment) {
   Int m = NTL::power(b, e) + c;
   m_b = b;
   m_e = e;
   m_c = c;
   init(m, decomp, filename, increment);
}

//===========================================================================

template<typename Int>
LCGComponent<Int>::~LCGComponent() {
   //a.kill();
   //orbitSeed.kill();
}

//============================================================================
template<typename Int>
void LCGComponent<Int>::init(const Int &m, DecompType decomp, const char *filename, bool increment) {
   m_m = m;
   m_increment = increment;
   setModulusIntP<Int>(m_m);
   if (increment) m_fact = IntFactorization<Int>(m_m);
   else m_fact = IntFactorization<Int>(m_m-1);
   m_fact.decompToFactorsInv (decomp, filename);
}

//===========================================================================

template<typename Int>
bool LCGComponent<Int>::maxPeriod(const Int &a) {
   if (!m_increment)
      return isPrimitiveElement(a, m_fact, m_m);   // Here, m is assumed to be prime.
   else {
      auto list = m_fact.getFactorList();
      for (auto iter = list.begin(); iter != list.end(); iter++) {
         if ((a - Int(1)) % (*iter).getFactor() != 0) return false;
      }
      if (m_m % 4 == 0 && (a - Int(1)) % 4 != 0) return false;
      return true;
   }
}

//===========================================================================

template<typename Int>
std::string LCGComponent<Int>::toString() {
   std::ostringstream os;
   os << "LCGComponent with: ";
   os << "\n   m = " << m_m << " = " << m_b << "^" << m_e << " + " << m_c;
   os << "\n   increment c_0 = " << m_increment << "\n";
   std::string str(os.str());
   return str;
}

template class LCGComponent<std::int64_t> ;
template class LCGComponent<NTL::ZZ> ;

} // End namespace LatMRG

#endif

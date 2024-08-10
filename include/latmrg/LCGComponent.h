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
 * The object is constructed for a given value of \f$m\f$, and can be used for several values of the
 * parameters \f$a\f$ and  \f$c_0\f$.  It does not look at the lattice structure.
 */
template<typename Int>
class LCGComponent {

//  using namespace LatMRG;

//private:
   // typedef NTL::vector<Int> IntVec;
   // typedef NTL::matrix<Int> IntMat;

public:

   /**
    * Constructor with modulus \f$m\f$, with the prime factor decomposition
    * of \f$m-1\f$ specified by `decompm1`.
    * The type `DecompType` is defined in `EnumTypes` and offers the
    * following choices:  `DECOMP`: the integer will be factored,
    * `DECOMP_WRITE`: it will be factored and the prime factors written in a file,
    * `DECOMP_READ`: the prime factors are read from a file,
    * `DECOMP_PRIME`: the integer is assumed to be prime.
    * The name `filem1` specifies the file where the factors of
    * \f$m-1\f$ are read or written when the type
    * <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt> is used.
    * The given files must be accessible by the program.
    *
    * Maybe we should just pass an `IntFactorization` instead of this.  Would be simpler?
    */
   LCGComponent(const Int &m, DecompType decompm1, const char *filem1);

   /**
    * Same as the previous constructor, but with `m=b^e+c`.
    */
   LCGComponent(const Int &b, int e, const Int &c, DecompType decompm1, const char *filem1);

   /**
    * For an LCG component with an additive term \f$c_0\f$, we also need the prime factor
    * decomposition of `m`.

   LCGComponent(const Int &m, DecompType decompm1, const char *filem1,
         DecompType decompm = NO_DECOMP, const char *filem = NULL);

   LCGComponent(const Int &b, int e, const Int &c, DecompType decompm1, const char *filem1,
         DecompType decompm = NO_DECOMP, const char *filem = NULL);
    */

   /**
    * Constructor similar to the above, except that the modulus of
    * congruence \f$m\f$ is inside the object `modul`.
    */
   // LCGComponent(Modulus<Int> &modul, int k, DecompType decom1,
   //        const char *filem1, DecompType decor, const char *filer);
   /**
    * Destructor.
    */
   ~LCGComponent();

   /**
    * Copy constructor;
    */
   // LCGComponent(const MRGComponent<Int> &comp);
   /**
    * Assignment operator.
    */
   // LCGComponent<Int>& operator=(const MRGComponent<Int> &comp);
   /**
    * Sets the multiplier \f$a\f$ of the recurrence to `a`.
    */
   void seta(const Int &a);

   /**
    * Gets the multiplier \f$a\f$ of the recurrence.
    */
   Int geta() const {
      return m_a;
   }

   /**
    * Returns `true` iff the LCG defined by
    * \f[
    *   x_n = a x_{n-1} \bmod m
    * \f]
    * has full period \f$m\f$.
    * This occurs if and only if \f$m\f$ is prime and \f$a\f$
    * is a primitive element modulo \f$m\f$.
    */
   bool maxPeriod(const Int &a);

   /**
    * Returns `true` iff the LCG defined by
    * \f[
    *   x_n = (a x_{n-1} + c_0) \bmod m
    * \f]
    * has full period when \f$c_0\f$ is relatively prime with \f$m\f$.
    * This function checks the following two conditions:
    * (1) Every prime divisor \f$q\f$ of \f$m\f$ must divide \f$a-1\f$;
    * (2) If 4 divides \f$m\f$, then it must divide \f$a-1\f$.
    */
   bool maxPeriodAdd(const Int &a);

   /**
    * Returns the value of the modulus \f$m\f$ of the recurrence.
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
    * Returns a const reference to the initial state `orbitSeed`, when available.
    * */
   Int& getOrbitSeed() {
      return orbitSeed;
   }

   /**
    * Returns this object as a string.
    */
   std::string toString();

   /**
    * The modulo of the MWC generator if we study one. This is because the
    * we check MWC period with module.m as the modulo of the equivalent LCG.
    * */
   Int m_MWCb;

private:

   /**
    * This is called by the constructor, with the same arguments.
    */
   void init(const Int &m, DecompType decompm1, const char *filem1);

   /**
    * Used for the prime factor decomposition of \f$m-1\f$.
    */
   IntFactorization<Int> ifm1;

   /**
    * The prime factor decomposition of \f$m\f$.
    * This one is used only for LCGs with a constant term.
    * Also used to compute the full period of a LCG with a carry.  ?????
    */
   IntFactorization<Int> ifm;

   /**
    * The modulus \f$m\f$ of the recurrence.
    */
   Int m_m;

   /**
    * The multiplier \f$a\f$ of the recurrence.
    */
   Int m_a;

   /**
    * Basis `b` when `m = b^e + c`.
    * */
   Int m_b;

   /**
    * Exponent `e` when `m = b^e + c`.
    * */
   int m_e;

   /**
    * Rest `c` when `m = b^e + c`.
    * */
   Int m_c;

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
   Int orbitSeed;

};   // End of class declaration


//===========================================================================
// IMPLEMENTATION

/*
 // Copy constructor.
 template<typename Int>
 LCGComponent<Int>::LCGComponent(const LCGComponent<Int> &lat) :
 ifm1(lat.ifm1), ifr(lat.ifr), m_k(lat.m_k) {
 module = lat.module;
 nj = lat.nj;
 rho = lat.rho;
 //   a.kill();
 m_a.resize(m_k);
 m_a = lat.m_a;
 //   orbitSeed.kill();
 orbitSeed.resize(m_k);
 orbitSeed = lat.orbitSeed;

 m_type = lat.m_type;
 if (m_type == MWC)
 m_MWCb = lat.m_MWCb;
 m_b = lat.m_b;
 m_e = lat.m_e;
 m_r = lat.m_r;
 }

 //===========================================================================

 // Assignment.
 template<typename Int>
 LCGComponent<Int>& LCGComponent<Int>::operator=(const LCGComponent<Int> &lat) {
 if (this != &lat) {
 m_k = lat.m_k;
 module = lat.module;
 nj = lat.nj;
 rho = lat.rho;
 //    a.kill();
 m_a.resize(m_k);
 LatticeTester::CopyVect(m_a, lat.m_a, m_k);
 //     orbitSeed.kill();
 orbitSeed.resize(m_k);
 LatticeTester::CopyVect(orbitSeed, lat.orbitSeed, m_k);
 //      ifm1 = lat.ifm1;
 //      ifr = lat.ifr;
 }
 m_type = lat.m_type;
 return *this;
 }
 */

//===========================================================================

template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &m, DecompType decompm1, const char *filem1) {
   m_b = m;
   m_e = 1;
   m_c = Int(0);
   init(m, decompm1, filem1);
}

//===========================================================================

template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &b, int e, const Int &c, DecompType decompm1, const char *filem1) {
   Int m = NTL::power(b, e) + c;
   m_b = b;
   m_e = e;
   m_c = c;
   init(m, decompm1, filem1);
}

//===========================================================================
/*
template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &m, DecompType decompm1, const char *filem1,
      DecompType decompm, const char *filem) {
   m_b = m;
   m_e = 1;
   m_c = Int(0);
   init(m, decompm1, filem1);
}

template<typename Int>
LCGComponent<Int>::LCGComponent(const Int &b, int e, const Int &c, DecompType decompm1, const char *filem1,
      DecompType decompm, const char *filem) {
   Int m = NTL::power(b, e) + c;
   m_b = b;
   m_e = e;
   m_c = c;
   init(m, decompm1, filem1);
}
*/

//===========================================================================

template<typename Int>
LCGComponent<Int>::~LCGComponent() {
   //a.kill();
   //orbitSeed.kill();
}

//============================================================================
template<typename Int>
void LCGComponent<Int>::init(const Int &m, DecompType decompm1, const char *filem1) {
   m_m = m;
   setModulusIntP<Int>(m_m);
   ifm = IntFactorization<Int>(m_m);
   // std::cout << "Inside init 3 \n";
   // ifm.factorizePlus();   // We always factorize m.
   ifm1 = IntFactorization<Int>(m_m-1);
   // ifm1(m_m - 1);
   if (decompm1 != NO_DECOMP) {
      if (decompm1 == DECOMP_READ) ifm1.read(filem1);
      else ifm1.factorizePlus();
      if (decompm1 == DECOMP_WRITE) {
         std::ofstream fout(filem1);
         fout << ifm1.toString();
      }
      // std::cout << "Inside init 5 \n";
      std::ofstream fout("dummy");
      // std::cout << "Before ifm1.toString \n";
      if (filem1 != 0) fout << ifm1.toString();
      // std::cout << "After ifm1.toString \n";
      remove("dummy");
   }
}

//===========================================================================

template<typename Int>
void LCGComponent<Int>::seta(const Int &a) {
   m_a = a;
}

//===========================================================================

template<typename Int>
bool LCGComponent<Int>::maxPeriod(const Int &a) {
   return isPrimitiveElement(a, ifm1, m_m);
}

//===========================================================================

template<typename Int>
bool LCGComponent<Int>::maxPeriodAdd(const Int &a) {
   auto list = ifm.getFactorList();
   for (auto iter = list.begin(); iter != list.end(); iter++) {
      if ((a - Int(1)) % (*iter).getFactor() != 0) return false;
   }
   if (m_m % 4 == 0 && (a - Int(1)) % 4 != 0) return false;
   return true;
}

//===========================================================================

template<typename Int>
std::string LCGComponent<Int>::toString() {
   std::ostringstream os;
   os << "LCGComponent:";
   os << "\n   m = " << m_m << " = " << m_b << "^" << m_e << " + " << m_c;
   os << "\n   a = " << m_a << "\n";
   std::string str(os.str());
   return str;
}

template class LCGComponent<std::int64_t> ;
template class LCGComponent<NTL::ZZ> ;

} // End namespace LatMRG

#endif

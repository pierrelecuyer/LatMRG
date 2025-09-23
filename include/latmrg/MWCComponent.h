#ifndef LATMRG_MWCCOMPONENT_H
#define LATMRG_MWCCOMPONENT_H

// #include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/LCGComponent.h"
#include <string>


namespace LatMRG {

using namespace LatMRG;

/**
 * This class represents a Multiply-with-carry (MWC) random number generator,
 * defined as in the guide by the recurrence:
 * \f{align}{
 *    \tau_n & = (a_1 x_{n-1} + \cdots + a_k x_{n-k} + c_{n-1}), \\
 *    x_n & = a_0^* \tau_n \bmod b, \\
 *    c_n & = \lfloor (\tau_n - a_0 x_n) / b \rfloor, \\
 * \f}
 * where gcd\f$(a_0, b)=1\f$, with the output
 * \f{align}{
 *    u_n & = \sum_{\ell=1}^\infty x_{n+\ell-1} b^{-\ell}.
 * \f}
 * In practice, the latter sum is truncated to just a few terms,
 * often just a single term.
 * This MCW generator has been shown to be equivalent to an LCG
 * with modulus \f$m = -a_0 + \sum_{j=1}^k a_j b^j $\f and multiplier
 * \f$b^*\f$, where \f$b^*\f$ is the multiplicative inverse of \f$b\f$
 * modulo \f$m\f$. This is the same as an LCG with multiplier \f $b\f$
 * with the outputs generated in reverse order.  The period length and
 * lattice structure are the same for both.
 * The aim of the present class is mainly to compute the parameters of the
 * corresponding LCG and also the period length of the MWC generator.
 * One can then analyze the lattice structure via `LCGLattice`.
 */
template<typename Int>
class MWCComponent {

public:

   /**
    * Constructs a `MWCComponent` with order `k`, modulus `b`, and vector of coefficients
    * `aa` of size `k+1` for which `a[j]` contains \f$a_j\f$ for \f$j=0,\dots,k\f$.
    * When not given, the modulus `b` and/or the vector `aa` must be set later via
    * `setb` or setbPow2, and `setaa`.  The order `k` cannot be changed.
    */
   MWCComponent(const int64_t k, const Int &b, const IntVec &aa);

   MWCComponent(const int64_t k, const Int &b);

   MWCComponent(const int64_t k);

   /**
    * Destructor.
    */
   ~MWCComponent();

   /**
    * Cleans and releases memory used by this object.
    */
   void kill();

   /**
    * Returns the order `k` of the recurrence.
    */
   int64_t getOrder() const {
      return m_k;
   }

   /**
    * Sets the vector of multipliers to `aa`, with `aa[j]` containing \f$a_j\f$.
    * The length of this vector must be the order \f$k\f$ of the MWC plus 1.
    */
   void setaa(const IntVec &aa);

   /**
    * Returns the vector of multipliers of the recurrence, in same format as `aa` above.
    */
   IntVec getaa() const {
      return m_aa;
   }

   /**
    * Sets the modulus `b`.
    */
   void setb(const Int &b);

   /**
    * Sets the modulus `b` to \f$b = 2^e\f$ and retain that it is a power of 2.
    */
   void setbPow2(const int64_t e);

   /**
    * Returns `true` iff `b` was set as a power of 2.
    */
   bool bIsPow2() { return (m_e > 0); }

   /**
    * Returns the modulus `b`.
    */
   Int getb() const {
      return m_b;
   }

   /**
    * Returns the inverse \f$b^*\f$, if it was computed before.
    */
   Int getLCGa() const {
      return m_LCGa;
   }

   /**
    * Returns the modulus \f$m\f$ of the corresponding LCG, if it was computed before.
    */
   Int getLCGm() const {
      return m_LCGm;
   }

   /**
    * Computes and returns the inverse \f$b^*\f$.
    */
   Int computeInvb();

   /**
    * Computes and returns the modulus \f$m\f$ of the corresponding LCG.
    */
   Int computeLCGModulus();

   /**
    * Returns the prime type of \f$m\f.
    * `numtrials` is the number of Miller-Rabin trials for probabilistic
    * primality testing.
    */
   PrimeType mIsPrime(const int64_t numtrials = 100);

   /**
    * Check if both \f$m\f and \f$(m-1)/2\f are prime.
    * This implies a maximal period of \f$(m-1)/2\f$ in case where
    * \f$a_0=1\f$ and \f$b\f$ is a power of 2,
    * `numtrials` is the number of Miller-Rabin trials for primality testing.
    */
   PrimeType mIsSafePrimePow2(const int64_t numtrials = 100);

   /**
    * Returns `true` if \f$m\f is prime or probably prime and if \f$b\f$
    * is a primitive root modulo \f$m\f$ (so the period is \f$m-1\f$.
    * `numtrials` is the number of Miller-Rabin trials for primality testing.
    */
   bool maxPeriod(const int64_t numtrials = 100);

   /**
    * Check if we have the maximal period of \f$(m-1)/2\f$ in case where \f$a_0=1\f$ and
    * \f$b\f$ is a power of 2.  This holds if \f$m\f is prime and 2 is a primitive
    * root modulo  \f$m\f$.
    * `numtrials` is the number of Miller-Rabin trials for primality testing.
    */
   bool maxPeriod1Pow2(const int64_t numtrials = 100);


private:

   /**
    * The order `k` of the MWC generator.
    */
   int m_k;

   /**
    * The modulo of the MWC generator this object represents.
    */
   Int m_b;

   /**
    * Used when `b = 2^e` is a power of 2.  Otherwise, `m_e = 0`.
    */
   int64_t m_e = 0;

   /**
    * The powers of `b`, from \f$b^0=1\f$ to \f$b^k\f$, precomputed
    * because we use them often.
    */
   IntVec m_powb;

   /**
    * The vector of coefficients \f$a_j\f$ of the MWC generator.
    * We should have `m_aa[j]`\f$a_j\f$, for \f$j\ge 0\f$.
    */
   IntVec m_aa;

   /**
    * The modulus of the corresponding LCG.
    */
   Int m_LCGm;

   /**
    * The multiplier of the corresponding LCG, which is `a = b^*`.
    */
   Int m_LCGa;

};
// End class declaration

//===========================================================================
// IMPLEMENTTION

//=============================================================================
// Main constructor.

template<typename Int>
MWCComponent<Int>::MWCComponent(const int64_t k, const Int &b, const IntVec &aa) :
      MWCComponent<Int>(k, b) {
   setaa(aa);
}

template<typename Int>
MWCComponent<Int>::MWCComponent(const int64_t k, const Int &b) :
      MWCComponent<Int>(k) {
   setb(b);
}

template<typename Int>
MWCComponent<Int>::MWCComponent(const int64_t k) {
   m_k = k;
   m_aa.SetLength(k + 1);
   m_powb.SetLength(k + 1);
}

//============================================================================

template<typename Int>
MWCComponent<Int>::~MWCComponent() {
   kill();
}

//===========================================================================

template<typename Int>
void MWCComponent<Int>::kill() {
   m_aa.kill();
   m_powb.kill();
}

//===========================================================================

template<typename Int>
void MWCComponent<Int>::setaa(const IntVec &aa) {
   assert(aa.length() == m_k + 1);
   assert(NTL::GCD(aa[0], m_b) == 1);
   m_aa = aa;
}

//===========================================================================

template<typename Int>
void MWCComponent<Int>::setb(const Int &b) {
   m_b = b;
   m_powb[0] = Int(1);
   for (int64_t j = 1; j <= m_k; j++) {
      m_powb[j] = m_powb[j-1] * b;
   }
}

//===========================================================================

template<typename Int>
void MWCComponent<Int>::setbPow2(const int64_t e) {
   m_e = e;
   setb(NTL::power(Int(2), e));
}


//============================================================================

// Computes `b^*`, the inverse of `b` modulo `m`.
template<typename Int>
Int MWCComponent<Int>::computeInvb() {
   if (m_aa[0] == Int(1))
      m_LCGa = (m_LCGm + 1) / m_b;
   else
      m_LCGa = NTL::InvMod(m_b, m_LCGm);
   return m_LCGa;
}

//============================================================================

template<typename Int>
Int MWCComponent<Int>::computeLCGModulus() {
   m_LCGm = -Int(m_aa[0]);
   for (int64_t j = 1; j <= m_k; j++) {
      m_LCGm += m_aa[j] * m_powb[j];
      //std::cout << "m_aa[j] = " << m_aa[j] << "\n";
   }
   return m_LCGm;
}

//============================================================================

template<typename Int>
PrimeType MWCComponent<Int>::mIsPrime(const int64_t numtrials) {
   return IntFactor<Int>::isPrime(m_LCGm, numtrials);
}

template<typename Int>
PrimeType MWCComponent<Int>::mIsSafePrimePow2(const int64_t numtrials) {
   assert(m_e > 0);
   assert(m_aa[0] == 1);
   return IntFactor<Int>::isSafePrime(m_LCGm, numtrials);
}

template<typename Int>
bool MWCComponent<Int>::maxPeriod(const int64_t numtrials) {
   PrimeType ptype = IntFactor<Int>::isPrime(m_LCGm, numtrials);
   if (ptype > 1) return false;
   IntFactorization<Int> m_fact(m_LCGm-Int(1));
   ptype = m_fact.decompToFactorsInv (DECOMP);
   if (ptype <= 1) return isPrimitiveElement(m_b, m_fact, m_LCGm);
   else return false;
}

template<typename Int>
bool MWCComponent<Int>::maxPeriod1Pow2(const int64_t numtrials) {
   //std::cout << "maxPeriod1Pow2 start \n";
   PrimeType ptype = IntFactor<Int>::isPrime(m_LCGm, numtrials);
   // std::cout << "  ptype = " << ptype << "\n";
   if (ptype > 1) return false;
   IntFactorization<Int> m_fact(m_LCGm-Int(1));
   // std::cout << "  m_fact created " << "\n";
   ptype = m_fact.decompToFactorsInv (DECOMP, NULL);
   // std::cout << "  decomp done, num factors: " <<  (m_fact.getFactorList()).size() << "\n";
   if (ptype <= 1) return isPrimitiveElement(Int(2), m_fact, m_LCGm);
   else return false;
}

//============================================================================

template class MWCComponent<std::int64_t> ;
template class MWCComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif


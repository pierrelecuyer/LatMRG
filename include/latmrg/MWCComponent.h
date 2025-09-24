#ifndef LATMRG_MWCCOMPONENT_H
#define LATMRG_MWCCOMPONENT_H

// #include "latticetester/EnumTypes.h"
// #include "latticetester/IntLatticeExt.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
// #include "latmrg/LCGComponent.h"
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
 *
 * The functions here handle MWC generators of order `k`, modulus `b`, and vector of coefficients
 * `aa` of size `k+1` for which `a[j]` contains \f$a_j\f$ for \f$j=0,\dots,k\f$.
 */

//===========================================================================
// Declarations

/**
 * Returns \f$b = 2^e\f$.
 */
template<typename Int>
static Int getPow2(int64_t e);

/**
 * Computes and returns the inverse \f$b^*\f$ of \f$b\f$ modulo \f$m\f$.
 */
template<typename Int>
static Int computeInvb(const Int &m, const Int &b);

/**
 * Also computes the inverse \f$b^*\f$, assuming that \f$b\f$ divides \f$m\f$.
 */
template<typename Int>
static Int computeInvbDivm(const Int &m, const Int &b);

/**
 * Computes and returns the modulus \f$m\f$ of the corresponding LCG.
 */
template<typename Int>
static Int computeLCGModulusMWC(const Int &b, const IntVec &aa);

/**
 * Returns the prime type of \f$m\f.
 * `numtrials` is the number of Miller-Rabin trials for probabilistic
 * primality testing.
 */
template<typename Int>
static PrimeType mIsPrime(const Int &m, const int64_t numtrials = 100);

/**
 * Check if both \f$m\f and \f$(m-1)/2\f are prime.
 * This implies a maximal period of \f$(m-1)/2\f$ in case where
 * \f$a_0=1\f$ and \f$b\f$ is a power of 2,
 * `numtrials` is the number of Miller-Rabin trials for primality testing.
 */
template<typename Int>
static PrimeType mIsSafePrime(const Int &m, const int64_t numtrials = 100);

/**
 * Returns `true` if \f$m\f is prime or probably prime and if \f$b\f$
 * is a primitive root modulo \f$m\f$ (so the period is \f$m-1\f$.
 * `numtrials` is the number of Miller-Rabin trials for primality testing.
 */
template<typename Int>
static bool maxPeriodMWC(const Int &m, const Int &b, const int64_t numtrials = 100);

/**
 * Check if we have the maximal period of \f$(m-1)/2\f$ in case where \f$a_0=1\f$ and
 * \f$b\f$ is a power of 2.  This holds if \f$m\f is prime and 2 is a primitive
 * root modulo  \f$m\f$.
 * `numtrials` is the number of Miller-Rabin trials for primality testing.
 */
template<typename Int>
static bool maxPeriod1Pow2MWC(const Int &m, const int64_t numtrials = 100);

// End class declaration

//===========================================================================
// IMPLEMENTTION

//=============================================================================

template<typename Int>
static Int getPow2(int64_t e) {
   return NTL::power(Int(2), e);
}

//============================================================================

// Computes `b^*`, the inverse of `b` modulo `m`.
template<typename Int>
static Int computeInvb(const Int &m, const Int &b) {
   return NTL::InvMod(b, m);
}

// Same, when we know that `b` divides `m` (e.g., if `a_0=1`).
template<typename Int>
static Int computeInvbDivm(const Int &m, const Int &b) {
   return (m + 1) / b;
}
//============================================================================

template<typename Int>
static Int computeLCGModulusMWC(const Int &b, const IntVec &aa) {
   assert(NTL::GCD(aa[0], b) == 1);
   int64_t k = aa.length();
   Int m(-aa[0]);
   Int powb(b);
   for (int64_t j = 1; j <= k; j++) {
      m += aa[j] * powb;
      powb *= b;
   }
   return m;
}

//============================================================================

template<typename Int>
static PrimeType mIsPrime(const Int &m, const int64_t numtrials) {
   return IntFactor < Int > ::isPrime(m, numtrials);
}

template<typename Int>
static PrimeType mIsSafePrime(const Int &m, const int64_t numtrials) {
   return IntFactor < Int > ::isSafePrime(m, numtrials);
}

template<typename Int>
static bool maxPeriodMWC(const Int &m, const Int &b, const int64_t numtrials) {
   PrimeType ptype = IntFactor < Int > ::isPrime(m, numtrials);
   if (ptype > 1) return false;
   IntFactorization < Int > fact(m - 1);
   ptype = fact.decompToFactorsInv(DECOMP);
   if (ptype <= 1) return isPrimitiveElement(b, fact, m);
   else return false;
}

template<typename Int>
static bool maxPeriod1Pow2MWC(const Int &m, const int64_t numtrials) {
   //std::cout << "maxPeriod1Pow2 start \n";
   PrimeType ptype = IntFactor < Int > ::isPrime(m, numtrials);
   // std::cout << "  ptype = " << ptype << "\n";
   if (ptype > 1) return false;
   IntFactorization < Int > fact(m - 1);
   // std::cout << "  m_fact created " << "\n";
   ptype = fact.decompToFactorsInv(DECOMP, NULL);
   // std::cout << "  decomp done, num factors: " <<  (m_fact.getFactorList()).size() << "\n";
   if (ptype <= 1) return isPrimitiveElement(Int(2), fact, m);
   else return false;
}

//============================================================================

// template class MWCComponent<std::int64_t> ;
// template class MWCComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif


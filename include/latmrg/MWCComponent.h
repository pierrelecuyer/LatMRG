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
 * defined as in the guide, by the recurrence:
 * \f{align}{
 *    \tau_n & = (a_1 x_{n-1} + \cdots + a_k x_{n-k} + c_{n-1}), \\
 *    x_n & = -a_0^* \tau_n \bmod b, \\
 *    c_n & = \lfloor (\tau_n + a_0 x_n) / b \rfloor, \\
 * \f}
 * where gcd\f$(a_0, b)=1\f$, with the output
 * \f{align}{
 *    u_n & = \sum_{\ell=1}^\infty x_{n+\ell-1} b^{-\ell}.
 * \f}
 * In practice, the latter sum is truncated to just a few terms,
 * often just a single term, but when we take the infinite sum,
 * this MCW generator is exactly equivalent to an LCG
 * with modulus \f$m = \sum_{j=0}^k a_j b^j $\f and multiplier
 * \f$b^*\f$, where \f$b^*\f$ is the multiplicative inverse of \f$b\f$
 * modulo \f$m\f$. This is the same as an LCG with multiplier \f $b\f$
 * with the outputs generated in reverse order. The period length and
 * lattice structure are the same for both.
 * The aim of the present class is mainly to compute the parameters of the
 * corresponding LCG and also the period length of the MWC generator.
 * One can then analyze the lattice structure via `LCGLattice`.
 * There is also a separate class `MWCLattice` that builds the basis in a
 * different way than in `LCGLattice`.
 *
 * The functions here handle MWC generators of order `k`, modulus `b`, and vector of coefficients
 * `aa` of size `k+1` for which `a[j]` contains \f$a_j\f$ for \f$j=0,\dots,k\f$.
 * The parameter `numtrials` is the number of Miller-Rabin trials for primality testing.
 *
 */

//===========================================================================
// Declarations

/**
 * Computes and returns the inverse \f$b^*\f$ of \f$b\f$ modulo \f$m\f$.
 */
template<typename Int>
static Int computeInvb(const Int &m, const Int &b);

/**
 * Also returns the inverse \f$b^*\f$, assuming that \f$b\f$ divides \f$m+1\f$.
 */
template<typename Int>
static Int computeInvbDivm(const Int &m, const Int &b);

/**
 * Computes and returns \f$a_0^*\f$, the inverse of \f$a_0\f$ modulo \f$b\f$.
 * The returned inverse will have the same sign as \f$a_0\f$.
 */
template<typename Int>
static Int computeInva0(const Int &b, const Int &a0);

/**
 * Computes and returns the modulus \f$m\f$ of the corresponding LCG.
 */
template<typename Int>
static Int computeLCGModulusMWC(const Int &b, const IntVec &aa);

/**
 * Returns the prime type of \f$m\f. For large numbers, the function performs a
 * Miller-Rabin probabilistic primality testing with `numtrials` trials.
 */
template<typename Int>
static bool mIsPrime(const Int &m, const int64_t numtrials = 100);

/**
 * Check if both \f$m\f and \f$(m-1)/2\f are prime.
 * If true, the returned PrimeType is less than 2.
 * This implies a period of at least \f$(m-1)/2\f$ in all cases,
 * and exactly \f$(m-1)/2\f$ if \f$b\f$ is an even power of 2.
 */
template<typename Int>
static bool mIsSafePrime(const Int &m, const int64_t numtrials = 100);

/**
 * Returns `true` if \f$m\f is prime or probably prime and if \f$b\f$
 * is a primitive root modulo \f$m\f$, so the period is \f$m-1\f$.
 * This is not possible if `b` is an even power of 2 or if \f$a_0 \bmod 8\f$
 * differs from 3 and 5. This function works for general values of `b` and `m`.
 */
template<typename Int>
static bool maxPeriodMWC(const Int &m, const Int &b, const int64_t numtrials = 100);

/**
 * Check if the period is at least \f$p = (m-1)/2\f$ when \f$b = 2^e\f$ and \f$m\f is prime.
 * This occurs if the order of 2 modulo \f$m\f is at least \f$p\f$ and gcd\f$(e, p) = 1\f$.
 * When `epow2 == true`, we assume that \f$e\f$ is a power of 2, so the latter holds iff \f$p\f$ is odd.
 * This function is useful for finding or testing cases for which \f$m\f is not a safe prime.
 */
template<typename Int>
static bool maxPeriodHalfMWC(const Int &m, const int64_t e, const bool epow2 = true, const int64_t numtrials = 100);

/**
 * For a MWC with parameters `b` and `aa`, and state \f$(x_{n-k+1},\dots,x_{n},c_n)\f$,
 * given in the vector `(xx[0..k])`, this function returns in `y` the LCG state \f$y_{n}\f$.
 * The order `k` is determined by the length of `aa`.
 */
template<typename Int>
static void MWCtoLCGState(Int &y, const Int &b, const IntVec &aa, const IntVec &xx);

/**
 * Same as the previous function, but this one returns in `y` the LCG state \f$y_{n-k}\f$.
 */
template<typename Int>
static void MWCtoLCGStateLagk(Int &y, const Int &b, const IntVec &aa, const IntVec &xx);

/**
 * This is the inverse of the previous function.
 * For a MWC with parameters `b` and `aa`, and LCG state \f$y_{n-k}\f$,
 * returns in `xx` the MWC state \f$(x_{n-k+1},\dots,x_{n},c_n)\f$.
 * The parameter `a0inv` provides the inverse of `aa[0]` modulo `b`.
 * If `a0inv` is not given, it is computed by the function. Should we remove this latter option?  ***
 */
template<typename Int>
static void LCGtoMWCStateLagk(IntVec &xx, const Int &b, const IntVec &aa, const Int &a0inv, const Int &y);

template<typename Int>
static void LCGtoMWCStateLagk(IntVec &xx, const Int &b, const IntVec &aa, const Int &y);

/**
 * This function is similar to the previous ones, except that it takes the LCG state
 * \f$y_{n}\f$ as input to compute the MWC state \f$(x_{n-k+1},\dots,x_{n},c_n)\f$.
 * It requires `m` in addition to `aa`.  The algorihm is different.
 */
template<typename Int>
static void LCGtoMWCStateDirect(IntVec &xx, const Int &b, const IntVec &aa, const Int &m, const Int &y);

/**
 * Returns the multiplier `jumpMult` required to jump ahead by `jumpSize` steps with the LCG.
 * The multiplier for one step ahead is the inverse of `b` modulo `m`.
 */
template<typename Int>
static void getJumpAheadMult(Int &jumpMult, const Int &b, Int &m, Int &jumpSize);

/**
 * Returns the multiplier `jumpMult` required to jump *backwards* by `jumpSize` steps with the LCG.
 * The multiplier to make one step backwards is `b`.
 */
template<typename Int>
static void getJumpBackMult(Int &jumpMult, const Int &b, Int &m, Int &jumpSize);

/**
 * Returns in `xx2` the MWC state which is `jumpSize` steps ahead of the current state in `xx` and `c`.
 * for a MWC with parameters `b`, `aa`, `m`, where `a0inv` is the inverse of `aa[0]` modulo `b`.
 * The vectors `xx` and `xx2` must have size `k+1` and contain a state \f$(x_0,\dots,x_{k-1},c)\f$.
 * This will jump ahead (using the inverse of `b`) if `jumpMult` was obtained via `getJumpMult`,
 * and will jump backwards if `jumpMult` was obtained via `getJumpBack`.
 * The transformations between the MWC and LCG states are done via `MWCtoLCGStateLagk` and
 * `LCGtoMWCStateLagk`. The variable `y` will hold the LCG state `y_{-k}` used in the process.
 * It permits one to recover this state after the jump, if desired. Its initial value is not used.
 */
template<typename Int>
static void jumpMWC(IntVec &xx2, const Int &b, const IntVec &aa, const Int &m,
      const Int &a0inv, IntVec &xx, Int &jumpMult, Int &y);

/**
 * Similar to `jumpMWC`, except that this function uses `MWCtoLCG` and `LCGtoMWCStateDirect`
 * for the transformations to/from the corresponding LCG state.
 */
template<typename Int>
static void jumpMWCDirect(IntVec &xx2, const Int &b, const IntVec &aa, const Int &m,
      IntVec &xx, Int &jumpMult, Int &y);

// End class declaration

//===========================================================================
// IMPLEMENTTION
//============================================================================

// Computes `b^*`, the inverse of `b` modulo `m`.
template<typename Int>
static Int computeInvb(const Int &m, const Int &b) {
   return NTL::InvMod(b, m);
}

// Special case for when we know that `b` divides `m+1` (e.g., if `a_0=-1`).
template<typename Int>
static Int computeInvbDivm(const Int &m, const Int &b) {
   return (m + 1) / b;
}

// Computes `a_0^*`, the inverse of `a_0` modulo `b`.
template<typename Int>
static Int computeInva0(const Int &b, const Int &a0) {
   Int a0inv;
   if (a0 == -1) a0inv = -1;
   else if (a0 == 1) a0inv = 1;
   else if (a0 < 0) a0inv = -NTL::InvMod(-a0, b);
   else a0inv = NTL::InvMod(a0, b);
   return a0inv;
}

template<typename Int>
static Int computeLCGModulusMWC(const Int &b, const IntVec &aa) {
   assert(NTL::GCD(aa[0], b) == 1);
   int64_t k = aa.length() - 1;
   Int m(aa[0]);
   // std::cout << "computeLCGModulusMWC, k = " << k << ",  m = " << m << "\n";
   Int powb(1);
   for (int64_t j = 1; j <= k; j++) {
      powb *= b;
      m += aa[j] * powb;
      // std::cout << " aa[j] = " << aa[j] << ", powb = " << powb << "\n";
   }
   return m;
}

//=========================================

template<typename Int>
inline static bool mIsPrime(const Int &m, const int64_t numtrials) {
   return NTL::ProbPrime(m, numtrials);
   // return IntFactor<Int>::isPrime(m, numtrials);
}

template<typename Int>
inline static bool mIsSafePrime(const Int &m, const int64_t numtrials) {
   if (!NTL::ProbPrime(m, numtrials)) return false;
   return NTL::ProbPrime((m - Int(1)) / Int(2), numtrials);
   // return IntFactor<Int>::isSafePrime(m, numtrials);
}

template<typename Int>
static bool maxPeriodMWC(const Int &m, const Int &b, const int64_t numtrials) {
   if (!IntFactor<Int>::isPrime(m, numtrials)) return false;  // m is not prime.
   IntFactorization<Int> fact(m - 1);
   PrimeType ptype = fact.decompToFactorsInv(DECOMP);
   if (ptype <= 1) return isPrimitiveElement(b, fact, m);
   else return false;
}

template<typename Int>
static bool maxPeriodHalfMWC(const Int &m, const int64_t e, const bool epow2, const int64_t numtrials) {
   if (!IntFactor<Int>::isPrime(m, numtrials)) return false;  // m is not prime.
   Int p = (m-1)/2;
   if ((epow2) && (p % 2) == 0) return false;  // epow2 and p is even
   if (NTL::GCD(conv<Int>(e), p) > 1) return false;
   IntFactorization<Int> fact(m-1);
   PrimeType ptype = fact.decompToFactorsInv(DECOMP, NULL);
   // std::cout << "  decomp done, num factors: " <  <  (m_fact.getFactorList()).size() << "\n";
   if (ptype <= 1) return isHighOrder(Int(2), p, fact, m);   // We only need order_m(2) >= p.
   else return false;
}

template<typename Int>
static void MWCtoLCGState(Int &y, const Int &b, const IntVec &aa, const IntVec &xx) {
   int64_t k = aa.length() - 1;
   Int bell = Int(1);  // b^\ell
   Int insum;   // inner sum
   y = xx[k];   // = c_n
   for (int64_t ell = 0; ell < k; ell++) {
      insum = Int(0);
      for (int64_t j = ell+1; j <= k; j++)
         insum += aa[j] * xx[ell-j+k];
      y += insum * bell;
      bell *= b;
   }
}

template<typename Int>
static void MWCtoLCGStateLagk(Int &y, const Int &b, const IntVec &aa, const IntVec &xx) {
   int64_t k = aa.length() - 1;
   Int bell = b;  // b^\ell
   Int insum;     // inner sum
   y = -aa[0] * xx[0];
   for (int64_t ell = 1; ell < k; ell++) {
      insum = Int(0);
      for (int64_t j = 0; j <= ell; j++)
         insum -= aa[j] * xx[ell-j];
      y += insum * bell;
      bell *= b;
   }
   y += xx[k] * bell;
}

template<typename Int>
static void LCGtoMWCStateLagk(IntVec &xx, const Int &b, const IntVec &aa, const Int &a0inv, const Int &y) {
   int64_t k = aa.length() - 1;
   Int sigma = y;
   for (int64_t i = 0; i < k; i++) {
      for (int64_t j = 1; j <= i; j++)
         sigma += aa[j] * xx[i-j];
      xx[i] = (-a0inv * sigma) % b;
      sigma = (sigma + aa[0] * xx[i]) / b;
   }
   xx[k] = sigma;
}

template<typename Int>
static void LCGtoMWCStateLagk(IntVec &xx, const Int &b, const IntVec &aa, const Int &y) {
   // Int a0inv = computeInva0(b, aa[0]);
   LCGtoMWCStateLagk(xx, b, aa, computeInva0(b, aa[0]), y);
}

template<typename Int>
static void LCGtoMWCStateDirect(IntVec &xx, const Int &b, const IntVec &aa, const Int &m, const Int &y) {
   int64_t k = aa.length() - 1;
   Int sigma = y;
   for (int64_t i = 0; i < k-1; i++) {
      xx[k-1-i] = (b * sigma) / m;
      sigma = b * sigma - m * xx[k-1-i];   // = (b * sigma) % m
   }
   xx[0] = (b * sigma) / m;
   xx[k] = 0;
   MWCtoLCGState(sigma, b, aa, xx);
   xx[k] = y - sigma;   // Here we recover c_n.
}

template<typename Int>
static void getJumpAheadMult(Int &jumpMult, const Int &b, Int &m, Int &jumpSize) {
   Int invb = computeInvb(m, b);
   NTL::PowerMod (jumpMult, invb, jumpSize, m);
}

template<typename Int>
static void getJumpBackMult(Int &jumpMult, const Int &b, Int &m, Int &jumpSize) {
   NTL::PowerMod (jumpMult, b, jumpSize, m);
}

template<typename Int>
static void jumpMWC(IntVec &xx2, const Int &b, const IntVec &aa, Int &m,
       Int &a0inv, IntVec &xx, Int &jumpMult, Int &y) {
   MWCtoLCGStateLagk(y, b, aa, xx);
   NTL::MulMod(y, jumpMult, y, m);
   LCGtoMWCStateLagk(xx2, b, aa, a0inv, y);
}

template<typename Int>
static void jumpMWCDirect(IntVec &xx2, const Int &b, const IntVec &aa, Int &m,
       IntVec &xx, Int &jumpMult, Int &y) {
   MWCtoLCGState(y, b, aa, xx);
   NTL::MulMod(y, jumpMult, y, m);
   LCGtoMWCStateDirect(xx2, b, aa, m, y);
}

//============================================================================

// template class MWCComponent<std::int64_t> ;
// template class MWCComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif


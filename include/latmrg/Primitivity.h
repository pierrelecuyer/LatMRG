#ifndef LATMRG_PRIMITIVITY_H
#define LATMRG_PRIMITIVITY_H

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZ_pXFactoring.h"
#include "NTL/ZZ_pE.h"
#include "NTL/lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
#include "NTL/lzz_pXFactoring.h"

#include "latticetester/Util.h"
#include "latmrg/FlexModInt.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"

namespace LatMRG {

/**
 * This file provides static functions to test for the primitivity of an integer
 * or of a polynomial in a finite field.
 *
 * Suppose that \f$a\f$, \f$e\f$, and \f$p\f$ are integers, \f$p\f$ is a
 * prime number, and \f$a\f$ is relatively prime with \f$m=p^e\f$.
 * The smallest positive integer \f$\lambda(m,a)\f$ for which
 * \f$a^{\lambda(m,a)}= 1 \ \bmod\ m \f$ is called the order of \f$a\f$ modulo
 * \f$m\f$. Any \f$a\f$ which has the maximum possible order for a given
 * \f$m\f$ is called a primitive element modulo \f$m\f$.
 * The function `isPrimitiveElement` tests for this property.

 * The function `isPrimitive` tests for the primitivity of the characteristic polynomial
 * of an MRG, with coefficients in \f$\mathbb Z_m\f$.  This is done by checking the
 * conditions proposed by Alanen and Knuth \cite mALA64a, and restated in
 * \cite sLEC23s and in the LatMRG user's guide.
 * The polynomial arithmetic is done using NTL.  A polynomial
 * \anchor REF__Primitivity_eq_poly1
 * \f[
 *   f(z) = c_0 + c_1 z^1 + \cdots + c_{k-1} z^{k-1} + c_k z^k    \tag{eq.poly1}
 * \f]
 * of degree \f$k\f$ and integer coefficients \f$c_j \in \mathbb Z_m\f$
 * is represented in NTL by the vector \f$(c_0,c_1,\dots,c_k)\f$ of its coefficients,
 * which is directly accessible.
 * The characteristic polynomial of an MRG, usually written as
 * \anchor REF__Primitivity_eq_poly2
 * \f[
 *   P(z) = z^k - a_1 z^{k-1} - \cdots- a_{k-1} z - a_k,   \tag{eq.poly2}
 * \f]
 * must be put in the form {@link REF__Primitivity_eq_poly1 (eq.poly1)} to use NTL.
 * It has \f$c_k=1, $c_{k-1}= -a_1\f, \cdots, c_1 = -a_{k-1}, c_0 = - a_k\f$.
 *
 * We recall that the polynomial \f$f(z)\f$ in {@link REF__Primitivity_eq_poly1 (eq.poly1)}
 * with \f$c_k=1\f$ is primitive modulo \f$m\f$ if and only if the three following conditions are satisfied:
 * \anchor REF__Primitivity_isprimi
 * <dl> <dt>None</dt>
 * <dd>
 * \f$[(-1)^{k} c_0]^{(m-1)/q} \bmod m \neq1\f$ for each prime
 * factor \f$q\f$ of \f$m - 1\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$z^r \bmod(f(z),m) =  (-1)^{k} c_0 \bmod m\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$z^{r/q} \bmod(f(z), m) \f$ has positive degree for each prime
 * factor \f$q\f$ of \f$r\f$, with \f$1<q< r\f$;
 * </dd>
 * </dl> where \f$r = (m^k - 1)/(m - 1)\f$.
 * Condition 1 is the same as saying that \f$(-1)^{k} c_0\f$ is a
 * primitive root of \f$m\f$. Condition 3 is automatically satisfied
 * when \f$r\f$ is prime.
 *
 * The types used for the polynomial coefficients in depend on the choice
 * of `Int` and are determined by the template class `FlexModInt.h`.
 * The function `isPrimitiveElement` does no use `FlexModInt`.
 */

//===========================================================================
// Declarations
/**
 * Returns `true` iff `a` is a primitive element modulo \f$p^e\f$.
 * The prime factor decomposition of \f$p-1\f$ must be given in `fac`,
 * and the list of inverse factors in `fac` must be up to date.
 */
template<typename Int>
static bool isPrimitiveElement(const Int &a, const IntFactorization<Int> &fac, const Int &p,
      long e = 1);

/**
 * Returns `true` iff the polynomial \f$f\f$ is a primitive polynomial
 * modulo \f$m\f$. The factorizations of \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
 * The modulus `m` should be set before, via the `setModulusIntP` function from `FlexModInt`.
 */
template<typename Int>
static bool isPrimitive(const NTL::Vec<Int> &aa, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr);

/**
 * Similar to `isPrimitive` above, except that this function only checks for Conditions 2 and 3.
 * Condition 1 is assumed to hold.
 */
template<typename Int>
static bool isPrimitive23(const NTL::Vec<Int> &aa, const Int &m,
      const IntFactorization<Int> &fr);

/**
 * Sets the coefficients of the polynomial `f` so it corresponds to the characteristic polynomial
 * \f$P(z) = z^k - a_1 z^{k-1} - \cdots- a_{k-1} z - a_k\f$ with coefficients \f$c_{k-j} = a_j = \f$`aa[j]`.
 * The modulus used in `IntP` is assumed to be correct, this is not verified.
 */
template<typename Int>
static void setCharacPoly(typename FlexModInt<Int>::PolX &f, const NTL::Vec<Int> &aa);

/**
 * In this version, the vector of coefficients is passed directly as an `IntVecP`,
 * so there is no need to convert it internally (this is more efficient).
 */
template<typename Int>
static void setCharacPoly(typename FlexModInt<Int>::PolX &f,
      typename FlexModInt<Int>::IntVecP &aaP);

/**
 * Converts the vector state `xx` of an MRG to its polynomial representation in `f`.
 * The vector of coefficients of the MRG is passed in `aa` and the modulus used in `IntP`
 * is assumed to be the correct one.
 * The vector `xx` is assumed to contain the vector state \f$(x_{n-k+1,\dots,x_n)\f$ with
 * `xx[j]` \f$x_{n-k+1+j]\f$.  The polynomial `f` will have degree `k-1` or less.
 */
template<typename Int>
static void vecMRGToPoly(typename FlexModInt<Int>::PolX &f, typename FlexModInt<Int>::IntVecP aaP,
      const NTL::Vec<Int> &xx);

/**
 * Converts the polynomial representation in `f` to the vector state `xx` of an MRG,
 * with `xx[j]` \f$x_{n-k+1+j]\f$.   The polynomial `f` must have degree `k-1` or less.
 * The vector of coefficients of the MRG is passed in `aa` and the modulus used in `IntP`
 * is assumed to be the correct one.
 */
template<typename Int>
static void polyToVecMRG(const NTL::Vec<Int> &xx, typename FlexModInt<Int>::IntVecP aaP,
      typename FlexModInt<Int>::PolX &f);

/**
 * Returns in `fj` the polynomial \f$x^j \mod f(x) (\bmod m)\f$ where \f$f(x)\f$ is the
 * polynomial with coefficients in 'C', modulo 'm'.
 */
//template<typename Int>
//static void powerMod(const Int &j, const IntVec &C, const Int &m, FlexModInt<Int>::PolX fj);

//===========================================================================
// IMPLEMENTTION

template<typename Int>
bool isPrimitiveElement(const Int &a, const IntFactorization<Int> &fac, const Int &p, long e) {
   if (0 == p) throw std::range_error("PrimitiveInt::isPrimitiveElement: p = 0");
   if (0 == a) return false;
   Int t1, t2;
   Int m;
   m = NTL::power(p, e);
   t1 = a;
   if (t1 < 0) t1 += m;
   const std::vector<Int> invList = fac.getInvFactorList();
   assert(!(invList.empty()));
   for (auto it = invList.begin(); it != invList.end(); it++) {
      if (*it == (m - 1)) continue;
      t2 = NTL::PowerMod(t1, *it, m);  // Works for either ZZ or int64_t.
      if (t2 == 1) return false;
   }
   return true;
}

template<typename Int>
static bool isPrimitive(const NTL::Vec<Int> &aa, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr) {
   int64_t k = aa.length() - 1;
   Int ak = aa[k];
   if ((k & 1) == 1) ak = -ak;
   if (!isPrimitiveElement(ak, fm, m)) return false;
   return isPrimitive23(aa, m, fr);
}

template<typename Int>
static bool isPrimitive23(const NTL::Vec<Int> &aa, const Int &m,
      const IntFactorization<Int> &fr) {
   // First, we make sure that IntP is using the correct modulus m.
   if (m != FlexModInt<Int>::IntP::modulus()) setModulusIntP<Int>(m);
   static typename FlexModInt<Int>::PolX f;
   setCharacPoly(f, aa);
   static int64_t k;
   k = deg(f);
   if (k == 1) return true;
   // First check irreducibility, which is faster.
   //  if (0 == DetIrredTest(Q))   // medium slow
   if (0 == IterIrredTest(f))   // fastest
   return false;

   // Test Condition 2
   std::cout << "Testing Condition 2 \n";
   Int r0;
   r0 = fr.getNumber();
   typename FlexModInt<Int>::PolX Q;
   Q = PowerXMod(r0, f);
   if (0 != deg(Q)) return false;
   Int T1;
   NTL::conv(T1, rep(ConstTerm(f)));
   if ((k & 1) == 1) T1 = -T1;
   if (T1 < 0) T1 += m;
   if (rep(ConstTerm(Q)) != T1) return false;

   // Test Condition 3
   std::cout << "Testing Condition 3 \n";
   if (fr.getStatus() == PRIME) return true;
   std::vector<Int> invFactorList = fr.getInvFactorList();
   assert(!invFactorList.empty());
   typename std::vector<Int>::const_iterator it = invFactorList.begin();
   while (it != invFactorList.end()) {
      Q = PowerXMod(*it, f);
      if (0 == deg(Q)) return false;
      ++it;
   }
   return true;
}

template<typename Int>
static void setCharacPoly(typename FlexModInt<Int>::PolX &f, const NTL::Vec<Int> &aa) {
   typename FlexModInt<Int>::IntVecP aaP;  // The coefficients `aa` must be converted to Z_p type.
   conv(aaP, aa);
   int64_t k = aaP.length() - 1;
   SetCoeff(f, k, 1);  // c_k = 1.
   for (int64_t j = 0; j < k; j++)
      SetCoeff(f, j, -aaP[k - j]);    // c_j = -a_{k-j} for j < k.
}

template<typename Int>
static void setCharacPoly(typename FlexModInt<Int>::PolX &f,
      typename FlexModInt<Int>::IntVecP aaP) {
   int64_t k = aaP.length() - 1;
   SetCoeff(f, k, 1);  // c_k = 1.
   for (int64_t j = 0; j < k; j++)
      SetCoeff(f, j, -aaP[k - j]);    // c_j = -a_{k-j} for j < k.
}

template<typename Int>
static void vecMRGToPoly(typename FlexModInt<Int>::PolX &f, typename FlexModInt<Int>::IntVecP aaP,
      const NTL::Vec<Int> &xx) {
   int64_t k = aaP.length() - 1;
   // We have xx[j] = x_{n-k+1+j].
   SetCoeff(f, k - 1, xx[0]);
   typename FlexModInt<Int>::IntP sum;
   for (int64_t j = 1; j < k; j++) {
      sum = xx[j];
      for (int64_t i = 1; i <= j; i++)
         sum -= aaP[i] * xx[j - i];
      SetCoeff(f, k - j, sum);
   }
}

template<typename Int>
static void polyToVecMRG(const NTL::Vec<Int> &xx, typename FlexModInt<Int>::IntVecP aaP,
      typename FlexModInt<Int>::PolX &f) {
   int64_t k = aaP.length() - 1;
   Getcoeff(xx[0], f, k - 1);
   typename FlexModInt<Int>::IntP sum;
   for (int64_t j = 1; j < k; j++) {
      Getcoeff(sum, f, k - j);
      for (int64_t i = 1; i <= j; i++)
         sum += aaP[i] * coeff(f, j - i);
      xx[j] = sum;
   }
}

//   power(fj, f, j);

}// namespace LatMRG

#endif

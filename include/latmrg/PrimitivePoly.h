#ifndef LATMRG_PRIMITIVEPOLY_H
#define LATMRG_PRIMITIVEPOLY_H

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
#include "latmrg/PrimitiveInt.h"

namespace LatMRG {

/*  Moved to FlexModInt.h
 *
template<typename Int>
class ModInt {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};

template<>
class ModInt<std::int64_t> {
public:
   typedef NTL::zz_p IntP;
   typedef NTL::zz_pE PolE;
   typedef NTL::vec_zz_p IntVecP;
   typedef NTL::mat_zz_p IntMatP;
   typedef NTL::zz_pX PolX;
};
*/

/*
template<>
class ModInt<NTL::ZZ> {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};
*/

/**
 * This file provides static functions to test the primitivity of polynomials
 * with coefficients in \f$\mathbb Z_m\f$.  This is done by checking the
 * conditions proposed by Alanen and Knuth \cite mALA64a, and restated in
 * \cite sLEC23s and in the LatMRG user's guide.
 * The polynomial arithmetic is done using NTL.  A polynomial
 * \anchor REF__PrimitivePoly_eq_poly1
 * \f[
 *   f(z) = c_0 + c_1 z^1 + \cdots + c_{k-1} z^{k-1} + c_k z^k    \tag{eq.poly1}
 * \f]
 * of degree \f$k\f$ and integer coefficients \f$c_j \in \mathbb Z_m\f$
 * is represented in NTL by the vector \f$(c_0,c_1,\dots,c_k)\f$ of its coefficients,
 * which is directly accessible.
 * The characteristic polynomial of an MRG, usually written as
 * \anchor REF__PrimitivePoly_eq_poly2
 * \f[
 *   P(z) = z^k - a_1 z^{k-1} - \cdots- a_{k-1} z - a_k,   \tag{eq.poly2}
 * \f]
 * must be put in the form {@link REF__PrimitivePoly_eq_poly1 (eq.poly1)} to use NTL.
 * It has \f$c_k=1, $c_{k-1}= -a_1\f, \cdots, c_1 = -a_{k-1}, c_0 = - a_k\f$.
 *
 * We recall that the polynomial \f$f(z)\f$ in {@link REF__PrimitivePoly_eq_poly1 (eq.poly1)}
 * with \f$c_k=1\f$ is primitive modulo \f$m\f$ if and only if the three following conditions are satisfied:
 * \anchor REF__PrimitivePoly_isprimi
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
 * The type `Int` is used to represent the polynomial coefficients.
 * The possible associated types `IntVec` are `int64_t*` and <tt>vec_ZZ</tt>.
 * The type `PolE` for the polynomials may be chosen
 * as <tt>zz_pE</tt> when \f$m\f$ is small enough, or may be set to
 * <tt>ZZ_pE</tt> which is implemented with the big integer type <tt>NTL::ZZ_p</tt>.
 */


/**
 * Sets to `m` the modulus used by NTL for its `IntP` calculations.
 */
//template<typename Int>
//static void setModulus(const Int &m);

/**
 * Returns `true` iff the polynomial \f$f\f$ is a primitive polynomial
 * modulo \f$m\f$. The factorizations of \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
 */
template<typename Int>
static bool isPrimitive(const NTL::vector<Int> &aa, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr);

/**
 * Similar to `isPrimitive` above, except that this function only checks for Conditions 2 and 3.
 * Condition 1 is assumed to hold.
 */
template<typename Int>
static bool isPrimitive23(const NTL::vector<Int> &aa, const Int &m,
      const IntFactorization<Int> &fr);

/**
 * Sets the coefficients of the polynomial `f` so it corresponds to the characteristic polynomial
 * \f$P(z) = z^k - a_1 z^{k-1} - \cdots- a_{k-1} z - a_k\f$ with coefficients \f$c_{k-j} = a_j = \f$`aa[j]`.
 */
template<typename Int>
static void setPoly(const NTL::vector<Int> &aa, typename ModInt<Int>::PolX &f);

/**
 * ****  This function does two different things !!!
 *   One is done by getPoly, one is done in NTL.   Remove.    ******
 * Returns in `fj` the polynomial \f$x^j \mod f(x) (\bmod m)\f$ where \f$f(x)\f$ is the
 * polynomial with coefficients in 'C', modulo 'm'.
 */
//template<typename Int>
//static void powerMod(const Int &j, const IntVec &C, const Int &m, ModInt<Int>::PolX fj);


//===========================================================================
// Implementation

/*
template<typename Int>
static void setModulus(const Int &m) {
   ModInt<Int>::IntP::init(m);        // Sets the modulus to m.
}
*/

template<typename Int>
static bool isPrimitive(const NTL::vector<Int> &aa, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr) {
   static typename ModInt<Int>::PolX f;  // We could also pass this f as a parameter and always reuse the same.
   setPoly(aa, f);
   static int64_t k;
   k = deg(f);
   Int c0;
   c0 = rep(ConstTerm(f));  // The constant term.
   if ((k & 1) == 1) c0 = -c0;
   if (!isPrimitiveElement(c0, fm, m)) return false;
   return isPrimitive23(aa, m, fr);
}

template<typename Int>
static bool isPrimitive23(const NTL::vector<Int> &aa, const Int &m,
      const IntFactorization<Int> &fr) {
   static typename ModInt<Int>::PolX f;
   setPoly(aa, f);
   static int64_t k;
   k = deg(f);
   if (k == 1) return true;
   // First test for irreducibility, which is faster.
   //  if (0 == DetIrredTest(Q))   // medium slow
   if (0 == IterIrredTest(f))   // fastest
   return false;

   // Test Condition 2
   std::cout << "Testing Condition 2 \n";
   Int r0;
   r0 = fr.getNumber();
   typename ModInt<Int>::PolX Q;
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
static void setPoly(const NTL::vector<Int> &aa, typename ModInt<Int>::PolX &f) {
   // ModInt<Int>::IntP::init(m);        // Sets the modulus to m.
   typename ModInt<Int>::IntVecP cP;  // The coefficients `aa` must be converted to Z_p type.
   conv(cP, aa);
   int64_t k = cP.length() - 1;
   SetCoeff(f, k, 1);
   for (int64_t i = 0; i < k; i++)
      SetCoeff(f, i, -cP[k - i]);
   f.normalize();  // Strip leading zeros.
}

//   power(fj, f, j);

}// namespace LatMRG
#endif

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
#include "latmrg/IntFactor.h"
#include "latmrg/PrimitiveInt.h"

namespace LatMRG {

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

template<>
class ModInt<NTL::ZZ> {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};

/**
 * This class deals with polynomials \f$P(x)\f$ in \f$\mathbb Z_m[X]\f$
 * defined as
 * \anchor REF__PrimitivePoly_eq_poly1
 * \f[
 *   P(x) = c_0 + c_1x^1 + c_2 x^2 + \cdots + c_n x^n \tag{eq.poly1}
 * \f]
 * with degree \f$n\f$ and integer coefficients \f$c_i\f$ in
 * \f$\mathbb Z_m\f$. The arithmetic operations on these polynomials are
 * done modulo \f$m\f$ and modulo a polynomial \f$f(x)\f$ of degree \f$k\f$.
 * Thus all polynomials will be reduced modulo \f$f(x)\f$. In LatMRG, the
 * modulus polynomial \f$f(x)\f$ is usually written in the form
 * \anchor REF__PrimitivePoly_eq_poly2
 * \f[
 *   f(x) = x^k - a_1x^{k-1} - \cdots- a_{k-1} x - a_k, \tag{eq.poly2}
 * \f]
 * and is associated with the recurrence
 * \anchor REF__PrimitivePoly_eq_rec2
 * \f[
 *   x_n = (a_1 x_{n-1} + a_2 x_{n-2} + \cdots+ a_k x_{n-k}) \bmod m. \tag{eq.rec2}
 * \f]
 * The two functions `setM` and `setF` *must* be called to initialize the
 * modulus \f$m\f$ and the modulus polynomial \f$f(x)\f$ before doing any
 * arithmetic operations on `PrimitivePoly` objects.
 *
 * We recall that the polynomial \f$f(x)\f$ in {@link REF__PrimitivePoly_eq_poly2 (eq.poly2)}
 * is a primitive polynomial modulo \f$m\f$ if and only if  the three
 * following conditions are satisfied:
 * \anchor REF__PrimitivePoly_isprimi
 * <dl> <dt>None</dt>
 * <dd>
 * \f$[(-1)^{k+1} a_k]^{(m-1)/q} \bmod m \neq1\f$ for each prime
 * factor \f$q\f$ of \f$m - 1\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$x^r \bmod(f(x),m) =  (-1)^{k+1} a_k \bmod m\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$x^{r/q} \bmod(f(x), m) \f$ has positive degree for each prime
 * factor \f$q\f$ of \f$r\f$, with \f$1<q< r\f$;
 * </dd>
 * </dl> where \f$r = (m^k - 1)/(m - 1)\f$.
 * Condition 1 is the same as saying that \f$(-1)^{k+1} a_k\f$ is a
 * primitive root of \f$m\f$. Condition 3 is automatically satisfied
 * when \f$r\f$ is prime.
 *
 * Type `Int` is used to represent the polynomial coefficients.
 * The possible associated types `IntVec` are `int64_t*`         // *****  Not NTL::vector<ZZ> ???
 *  and <tt>vec_ZZ</tt>. Type `PolE` for the polynomials may be chosen
 * as <tt>zz_pE</tt> when \f$m\f$ is small enough, or may be set to
 * <tt>ZZ_pE</tt> which is implemented with the big integer type <tt>NTL::ZZ_p</tt>.
 *
 */
template<typename Int>
class PrimitivePoly: public ModInt<Int>::PolE {

private:
   typedef NTL::vector<Int> IntVec;

public:

   /**
    * Minimal constructor: this object is set to the zero polynomial.
    */
   PrimitivePoly();

   /**
    * Initializes this object to the polynomial represented by \f$j\f$.   // ????
    */
   void setVal(int64_t j);

   /**
    * Initializes this object to the polynomial with coefficients in \f$C\f$.
    */
   void setVal(const IntVec &C);

   /**
    * Initializes this object to the polynomial represented in `str`.
    */
   void setVal(std::string &str);

   const typename ModInt<Int>::PolX& getVal() {
      return rep(*this);
   }

   /**
    * Sets this polynomial to \f$x^j \mod f(x) (\bmod m)\f$.
    */
   void powerMod(const Int &j);

   /**
    * Returns in `C` the coefficients of this polynomial as a vector of
    * \f$k\f$ components, where \f$k\f$ is the degree of the modulus \f$f(x)\f$.
    */
   void toVector(IntVec &C);

   /**
    * Initializes the modulus \f$m\f$ for this class. This must be called before
    * doing any operations on `PrimitivePoly` objects.
    */
   static void setM(const Int &m);

   /**
    * Returns a read-only reference to \f$m\f$.
    */
   static const Int& getM() {
      return m_m;
   }

   /**
    * Initializes the modulus polynomial
    * \f[
    *   f(x) = c_0 + c_1 x +c_2 x^2 + \cdots+ c_k x^k
    * \f]
    * of degree \f$k\f$ for this class from the
    * coefficients \f$c_i = \f$`C[i]` of vector `C` of dimension
    * \f$k + 1\f$. This function must be called before doing any
    * arithmetic operations on `PrimitivePoly` objects.
    */
   static void setF(const typename ModInt<Int>::IntVecP &C);

   /**
    * Same as above, but instead from an `IntVec` object.
    */
   static void setF(const IntVec &C);

   /**
    * Returns the modulus polynomial \f$f(x)\f$.
    */
   static const typename ModInt<Int>::PolX& getF() {
      return ModInt<Int>::PolE::modulus().val();
   }

   /**
    * Returns the degree of the modulus polynomial \f$f(x)\f$.
    */
   static int64_t getK() {
      return ModInt<Int>::PolE::degree();
   }

   /**
    * Given a vector \f$C = [c_0, c_1, \dots, c_{k-1}, c_k]\f$, this function
    * reorders the components of \f$C\f$ in the form \f$C = [c_k, c_{k-1},
    * ..., c_1, c_0]\f$ for `kind = 1`, and in the form \f$C =
    * [-c_k, -c_{k-1}, ..., -c_1, 1]\f$ for `kind = 2`. For other values of
    * `kind`, nothing is done.
    */
   static void reverse(IntVec &c, int64_t k, int64_t kind);

   /**
    * Returns `true` if the modulus \f$f(x)\f$ is a primitive polynomial
    * modulo \f$m\f$. The factorizations of
    * \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
    */
   bool isPrimitive(const PrimitiveInt<Int> &fm, const IntFactorization<Int> &fr);

   /**
    * This method returns `true` iff the primitivity conditions 2 and 3 given above
    * are satisfied by the modulus \f$f(x)\f$. It
    * does not check condition 1, assuming it to be `true`.
    */
   bool isPrimitive(const IntFactorization<Int> &r);

   /**
    * Returns a representation of this object as a string.
    */
   std::string toString() const;

private:

   /**
    * Modulus of the congruence.
    */
   static Int m_m;

   /**
    * Degree of the modulus polynomial \f$f\f$.
    */
   static int64_t m_k;

   /**
    * The polynomial \f$f(x) = x\f$.
    */
   static typename ModInt<Int>::PolX m_x;
};

template<typename Int>
Int PrimitivePoly<Int>::m_m;

template<typename Int>
int64_t PrimitivePoly<Int>::m_k;

template<typename Int>
typename ModInt<Int>::PolX PrimitivePoly<Int>::m_x;

/*=========================================================================*/
//PrimitivePoly::PrimitivePoly (const IntVec & C, int64_t n)
//{
//   init (C, n);
//}
/*=========================================================================*/

template<typename Int>
PrimitivePoly<Int>::PrimitivePoly() {
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setM(const Int &m) {
   ModInt<Int>::IntP::init(m);
   m_m = m;
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setF(const typename ModInt<Int>::IntVecP &c) {
   m_k = c.length() - 1;
   typename ModInt<Int>::PolX f;
   for (int64_t i = 0; i <= m_k; i++)
      SetCoeff(f, i, c[i]);
   f.normalize();
   ModInt<Int>::PolE::init(f);
   std::string str = "[0 1]";
   std::istringstream in(str);
   m_x.kill();
   in >> m_x;
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setF(const IntVec &V) {
   typename ModInt<Int>::IntVecP vP;
   conv(vP, V);
   setF(vP);
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setVal(int64_t j) {
   typename ModInt<Int>::IntP coeff;
   conv(coeff, 1);
   typename ModInt<Int>::PolX val;
   SetCoeff(val, j, coeff);
   val.normalize();
   ModInt<Int>::PolE::init(val);
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setVal(const IntVec &C) {
   std::ostringstream out;
   out << "[";
   for (int64_t i = 0; i < C.length(); i++) {
      out << C[i] << " ";
   }
   out << "]";
   std::string str = out.str();
   setVal(str);
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::powerMod(const Int &j) {
   std::string str = "[0 1]";
   std::istringstream in(str);
   PrimitivePoly<Int> A;
   in >> A;
   power(*this, A, j);
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::setVal(std::string &str) {
   std::istringstream in(str);
   in >> *this;
}

//===========================================================================

template<typename Int>
std::string PrimitivePoly<Int>::toString() const {
   std::ostringstream sortie;
   typename ModInt<Int>::PolX pX = ModInt<Int>::PolE::modulus().val();
   sortie << "m = " << PrimitivePoly<Int>::m_m << std::endl;
   sortie << "k = " << getK() << std::endl;
   sortie << "f = " << getF() << std::endl;
   sortie << "v = " << *this << std::endl;
   sortie << "    " << std::endl;
   return sortie.str();
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::toVector(IntVec &C) {
   for (int64_t i = 0; i < getK(); ++i) {
      C[i] = rep(coeff(rep(*this), i));
   }
}

/*=========================================================================*/
#if 0
  bool PrimitivePoly<Int>::isIrreducible ()
  {
    // Implements the Crandall-Pomerance method.  Its speed is nearly identical
    // to that of the `NTL::DetIrredTest` method, and slower than `NTL::IterIrredTest`
    // especially for large `k`.  So we use the latter.
    PolE g;
    PolX d;
    std::string str = "[0 1]";
    std::istringstream in (str);
    in >> g;
    for (int64_t i = 1; i <= m_k/2; ++i) {
      power(g, g, m_m);
      GCD (d, PolE::modulus(), rep(g) - m_x);
      if (1 != IsOne(d))
        return false;
    }
    return true;
  }
#endif

/*=========================================================================*/

template<typename Int>
bool PrimitivePoly<Int>::isPrimitive(const IntFactorization<Int> &r) {
   if (1 == getK()) return true;
   // Is f irreducible ?
   typename ModInt<Int>::PolX Q;
   Q = getF();
   //  if (!isIrreducible())      // slow
   //  if (0 == DetIrredTest(Q))   // medium slow
   if (0 == IterIrredTest(Q))   // fastest
   return false;

   // ---- Test Condition 2
   Int r0;
   r0 = r.getNumber();
   Q = PowerXMod(r0, getF());
   if (0 != deg(Q)) return false;
   Int T1;
   T1 = rep(ConstTerm(getF()));
   if ((getK() & 1) == 1) T1 = -T1;
   if (T1 < 0) T1 += getM();
   if (rep(ConstTerm(Q)) != T1) return false;

   // ---- Test Condition 3
   if (r.getStatus() == LatticeTester::PRIME) return true;
   std::vector<Int> invFactorList = r.getInvFactorList();
   //   assert (!invFactorList.empty ());
   typename std::vector<Int>::const_iterator it = invFactorList.begin();

   while (it != invFactorList.end()) {
      Q = PowerXMod(*it, getF());
      if (0 == deg(Q)) return false;
      ++it;
   }
   return true;
}

/*=========================================================================*/

template<typename Int>
bool PrimitivePoly<Int>::isPrimitive(const PrimitiveInt<Int> &fm, const IntFactorization<Int> &fr) {
   Int a0;
   // rep is the NTL::ZZ equivalent of the NTL::ZZ_p element.
   a0 = -rep(ConstTerm(getF()));
   if ((getK() & 1) == 0) a0 = -a0;
   if (!fm.isPrimitiveElement(a0)) {
      return false;
   }
   return isPrimitive(fr);
}

/*=========================================================================*/

template<typename Int>
void PrimitivePoly<Int>::reverse(IntVec &C, int64_t n, int64_t kind) {
   int64_t i;
   Int temp;
   if ((kind == 2) || (kind == 1)) {
      for (i = 0; i < (n + 1) / 2; i++) {
         temp = C[n - i];
         C[n - i] = C[i];
         C[i] = temp;
      }
   }
   if (kind == 2) {
      for (i = 0; i < n; i++)
         C[i] = -C[i];
      C[n] = 1;
   }
}

template class PrimitivePoly<std::int64_t> ;
template class PrimitivePoly<NTL::ZZ> ;

}
#endif

#ifndef LATMRG_POLY_H
#define LATMRG_POLY_H

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
#include "latticetester/IntFactor.h"

#include "latmrg/IntPrimitivity.h"

namespace LatMRG {
  template<typename Int>
    class ModInt{
      public:
        typedef NTL::ZZ_p IntP;
        typedef NTL::ZZ_pE PolE;
        typedef NTL::vec_ZZ_p IntVecP;
        typedef NTL::mat_ZZ_p IntMatP;
        typedef NTL::ZZ_pX PolX;
    };
  template<> class ModInt<std::int64_t> {
    public:
      typedef NTL::zz_p IntP;
      typedef NTL::zz_pE PolE;
      typedef NTL::vec_zz_p IntVecP;
      typedef NTL::mat_zz_p IntMatP;
      typedef NTL::zz_pX PolX;
  };
  template<> class ModInt<NTL::ZZ> {
    public:
      typedef NTL::ZZ_p IntP;
      typedef NTL::ZZ_pE PolE;
      typedef NTL::vec_ZZ_p IntVecP;
      typedef NTL::mat_ZZ_p IntMatP;
      typedef NTL::ZZ_pX PolX;
  };

  /**
   * This class implements polynomials \f$P(x)\f$ in \f$\mathbb Z_m[X]\f$
   * defined as
   * \anchor REF__PolyPE_eq_poly1
   * \f[
   *   P(x) = c_0 + c_1x^1 + c_2 x^2 + \cdots+ c_nx^n \tag{eq.poly1}
   * \f]
   * with degree \f$n\f$ and integer coefficients \f$c_i\f$ in
   * \f$\mathbb Z_m\f$. The arithmetic operations on objects of this class are
   * done modulo \f$m\f$ and modulo a polynomial \f$f(x)\f$ of degree \f$k\f$.
   * Thus all polynomials will be reduced modulo \f$f(x)\f$. In LatMRG, the
   * modulus polynomial \f$f(x)\f$ is usually written in the form
   * \anchor REF__PolyPE_eq_poly2
   * \f[
   *   f(x) = x^k - a_1x^{k-1} - \cdots- a_{k-1} x - a_k, \tag{eq.poly2}
   * \f]
   * and is associated with the recurrence
   * \anchor REF__PolyPE_eq_rec2
   * \f[
   *   x_n = (a_1x_{n-1} + a_2x_{n-2} + \cdots+ a_k x_{n-k}) \bmod m. \tag{eq.rec2}
   * \f]
   * The two functions `setM` and `setF` *must* be called to initialize the
   * modulus \f$m\f$ and the modulus polynomial \f$f(x)\f$ before doing any
   * arithmetic operations on `PolyPE` objects, otherwise the results are
   * unpredictable.
   *
   * Type `Int` is used to represent polynomial coefficients. It may be
   * chosen as `long` for \f$m < 2^{50}\f$ (on 64-bit machines), or as the big
   * integer type `ZZ` otherwise. The possible associated types `IntVec` are
   * `long*` and <tt>vec_ZZ</tt>. Type `PolE` for the polynomials may be chosen
   * as <tt>zz_pE</tt> when \f$m < 2^{50}\f$, or it may be set to
   * <tt>ZZ_pE</tt> which is implemented with the big integer type
   * <tt>ZZ_p</tt>.
   *
   */
  template<typename Int>
    class PolyPE : public ModInt<Int>::PolE {
      private:
        typedef NTL::vector<Int> IntVec;
      public:

        /**
         * Initializes the modulus \f$m\f$ for this class. This must be called before
         * doing any operations on `PolyPE` objects, otherwise the results are
         * unpredictable.
         */
        static void setM (const Int & m);

        /**
         * Returns a read-only reference to \f$m\f$.
         */
        static const Int & getM () { return m_m; }

        /**
         * Initializes the modulus polynomial \f$f(x) = c_0 + c_1x +c_2x^2 +
         * \cdots+ c_kx^k\f$ of degree \f$k\f$ for this class from the
         * coefficients \f$c_i = \f$<tt>C[i]</tt> of vector `C` of dimension
         * \f$k + 1\f$. This function must be called before doing any
         * arithmetic operations on `PolyPE` objects, otherwise the results are
         * unpredictable.
         */
        static void setF (const typename ModInt<Int>::IntVecP & C);

        /**
         * Same as above, but instead from a `IntVec`.
         */
        static void setF (const IntVec & C);

        /**
         * Returns the modulus polynomial \f$f(x)\f$.
         */
        static const typename ModInt<Int>::PolX & getF () { return ModInt<Int>::PolE::modulus().val(); }

        /**
         * Returns the degree of the modulus \f$f(x)\f$.
         */
        static long getK () { return ModInt<Int>::PolE::degree(); }

        /**
         * Given a vector \f$C = [c_0, c_1, …, c_{k-1}, c_k]\f$, this function
         * reorders the components of \f$C\f$ in the form \f$C = [c_k, c_{k-1},
         * …, c_1, c_0]\f$ for `kind = 1`, and in the form \f$C =
         * [-c_k, -c_{k-1}, …, -c_1, 1]\f$ for `kind = 2`. For other values of
         * `kind`, it has no effect.
         */
        static void reverse (IntVec & c, long k, int kind);

        /**
         * Minimal constructor: this object is set to the **0** polynomial.
         */
        PolyPE ();
        const typename ModInt<Int>::PolX & getVal () { return rep(*this); }
        void setVal (long j);

        /**
         * Initializes this object to \f$C\f$.
         */
        void setVal (const IntVec & C);

        /**
         * Initializes this object to the polynomial in `str`.
         */
        void setVal (std::string & str);

        /**
         * Sets \f$v = x^j \mod f(x) (\bmod m)\f$.
         */
        void powerMod (const Int & j);

        /**
         * Returns the coefficients of this polynomial as a vector \f$C\f$ of
         * \f$k\f$ components, where \f$k\f$ is the degree of the modulus
         * \f$f(x)\f$.
         */
        void toVector (IntVec & c);

        /**
         * Returns `true` if the modulus \f$f(x)\f$ is a primitive polynomial
         * modulo \f$m\f$. For this to be true, assuming that \f$f(x)\f$ has
         * the form {@link REF__PolyPE_eq_poly2 (eq.poly2)} above, the three
         * following conditions must be satisfied: \anchor REF__PolyPE_isprimi
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
         * </dl> where \f$r = (m^k - 1)/(m - 1)\f$. The factorizations of
         * \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
         * Condition 1 is the same as saying that \f$(-1)^{k+1} a_k\f$ is a
         * primitive root of \f$m\f$. Condition 3 is automatically satisfied
         * when \f$r\f$ is prime.
         */
        bool isPrimitive (const IntPrimitivity<Int> & fm, const IntFactorization<Int> & fr);

        /**
         * Given the factorization of \f$r\f$, this method returns `true` if
         * conditions 2 and 3 above are satisfied by the modulus \f$f(x)\f$. It
         * does not check condition 1, assuming it to be `true`.
         */
        bool isPrimitive (const IntFactorization<Int> & r);

        /**
         * Returns this object as a string.
         */
        std::string toString () const;
      private:

        /**
         * Modulus of congruence.
         */
        static Int m_m;

        /**
         * Degree of the modulus polynomial \f$f\f$.
         */
        static long m_k;

        /**
         * The polynomial \f$f(x) = x\f$.
         */
        static typename ModInt<Int>::PolX m_x;
    };

  template<typename Int>
    Int PolyPE<Int>::m_m;
  template<typename Int>
    long PolyPE<Int>::m_k;
  template<typename Int>
    typename ModInt<Int>::PolX PolyPE<Int>::m_x;



  /*=========================================================================*/
  //PolyPE::PolyPE (const IntVec & C, long n)
  //{
  //   init (C, n);
  //}


  /*=========================================================================*/

  template<typename Int>
    PolyPE<Int>::PolyPE ()
    {}


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setM (const Int & m)
    {
      ModInt<Int>::IntP::init (m);
      m_m = m;
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setF (const typename ModInt<Int>::IntVecP & c)
    {
      m_k = c.length () - 1;
      typename ModInt<Int>::PolX f;
      for (long i = 0; i <= m_k; i++)
        SetCoeff (f, i, c[i]);
      f.normalize ();
      ModInt<Int>::PolE::init (f);
      std::string str = "[0 1]";
      std::istringstream in (str);
      m_x.kill();
      in >> m_x;
    }

  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setF (const IntVec & V)
    {
      typename ModInt<Int>::IntVecP vP;
      conv (vP, V);
      setF (vP);
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setVal (long j)
    {
      typename ModInt<Int>::IntP coeff;
      conv (coeff, 1);
      typename ModInt<Int>::PolX val;
      SetCoeff (val, j, coeff);
      val.normalize ();
      ModInt<Int>::PolE::init (val);
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setVal (const IntVec & C)
    {
      std::ostringstream out;
      out << "[";
      for (int i = 0; i < C.length (); i++) {
        out << C[i] << " ";
      }
      out << "]";
      std::string str = out.str ();
      setVal (str);
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::powerMod (const Int & j)
    {
      std::string str = "[0 1]";
      std::istringstream in (str);
      PolyPE<Int> A;
      in >> A;
      power (*this, A, j);
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::setVal (std::string & str)
    {
      std::istringstream in (str);
      in >> *this;
    }


  //===========================================================================

  template<typename Int>
    std::string PolyPE<Int>::toString ()const
    {
      std::ostringstream sortie;

      typename ModInt<Int>::PolX pX = ModInt<Int>::PolE::modulus ().val ();
      sortie << "m = " << PolyPE<Int>::m_m << std::endl;
      sortie << "k = " << getK () << std::endl;
      sortie << "f = " << getF () << std::endl;
      sortie << "v = " << *this << std::endl;
      sortie << "    " << std::endl;
      return sortie.str ();
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::toVector (IntVec & C)
    {
      // C = rep(*this);

      // remark: in the code it was "i <= getK()" before but this produces
      // errors because sometimes dim(C) < getK().

      for (int i = 0; i < getK (); ++i) {
        C[i] = rep (coeff (rep (*this), i));
      }
    }


  /*=========================================================================*/
#if 0
  bool PolyPE<Int>::isIrreducible ()
  {
    // Méthode de crandall-Pomerance. La vitesse de cette méthode est
    // pratiquement identique à celle de NTL::DetIrredTest, mais plus
    // lente que NTL::IterIrredTest, spécialement pour grand k.
    PolE g;
    PolX d;
    std::string str = "[0 1]";
    std::istringstream in (str);
    in >> g;
    for (int i = 1; i <= m_k/2; ++i) {
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
    bool PolyPE<Int>::isPrimitive (const IntFactorization<Int> & r)
    {
      if (1 == getK())
        return true;

      // Is f irreducible ?
      typename ModInt<Int>::PolX Q;
      Q = getF ();
      //  if (!isIrreducible())      // slow
      //  if (0 == DetIrredTest(Q))   // medium slow
      if (0 == IterIrredTest(Q))   // fastest
        return false;

      // ---- Condition 2
      Int r0;
      r0 = r.getNumber ();
      Q = PowerXMod (r0, getF ());
      if (0 != deg (Q))
        return false;

      Int T1;
      T1 = rep (ConstTerm (getF ()));
      if ((getK () & 1) == 1)
        T1 = -T1;
      if (T1 < 0)
        T1 += getM ();
      if (rep (ConstTerm (Q)) != T1)
        return false;

      // ---- Condition 3
      if (r.getStatus () == LatticeTester::PRIME)
        return true;

      std::vector <Int> invFactorList = r.getInvFactorList ();
      //   assert (!invFactorList.empty ());
      typename std::vector <Int>::const_iterator it = invFactorList.begin ();

      while (it != invFactorList.end ()) {
        Q = PowerXMod (*it, getF());
        if (0 == deg (Q))
          return false;
        ++it;
      }
      return true;
    }


  /*=========================================================================*/

  template<typename Int>
    bool PolyPE<Int>::isPrimitive (const IntPrimitivity<Int> & fm,
        const IntFactorization<Int> & fr)
    {
      Int a0;
      a0 = -rep (ConstTerm (getF ()));
      if ((getK () & 1) == 0)
        a0 = -a0;
      if (!fm.isPrimitiveElement (a0)) {
        return false;
      }
      return isPrimitive (fr);
    }


  /*=========================================================================*/

  template<typename Int>
    void PolyPE<Int>::reverse (IntVec & C, long n, int kind)
    {
      long i;
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

  extern template class PolyPE<std::int64_t>;
  extern template class PolyPE<NTL::ZZ>;

}
#endif

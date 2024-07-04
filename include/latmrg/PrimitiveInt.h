#ifndef LATMRG_PRIMITIVEINT_H
#define LATMRG_PRIMITIVEINT_H

// #include "latticetester/ntlwrap.h"

#include "latmrg/IntFactorization.h"

namespace LatMRG {

  /**
   * This class deals with primitive roots and primitive elements modulo an
   * integer. Suppose that \f$a\f$, \f$e\f$, and \f$p\f$ are integers, \f$p\f$ is a
   * prime number, and \f$a\f$ is relatively prime with \f$m=p^e\f$.
   * The smallest positive integer \f$\lambda(m,a)\f$ for which
   * \f$a^{\lambda(m,a)}= 1 \ \bmod\ m \f$ is called the order of \f$a\f$ modulo
   * \f$m\f$. Any \f$a\f$ which has the maximum possible order for a given
   * \f$m\f$ is called a *primitive root* modulo \f$m\f$.
   * This class provides functions to test for this property.
   *
   * *****   This is incomplete.  Must be clarified.  And say which algo is used!  ******
   *
   * For the following
   * important cases, the value of the order for given \f$m\f$ is \cite rKNU98a :
   * \f{align*}{
   *  \lambda(2^e)
   *    &
   *   =
   *    2^{e-2}, \qquad e > 2
   *  \\
   *  \lambda(p^e)
   *    &
   *   =
   *    p^{e-1}(p - 1), \qquad p > 2.
   * \f}
   */
  template<typename Int>  class PrimitiveInt {

      private:
        typedef NTL::vector<Int> IntVec;

      public:

        PrimitiveInt ();

        /**
         * Constructor that fixes the modulus of congruence to \f$m = p^e\f$.
         * The argument `f` must contain the prime factor decomposition of
         * \f$p-1\f$ with its inverse factors also precomputed.
         */
        PrimitiveInt (const IntFactorization<Int> & f, const Int & p, long e = 1);

        /**
         * Returns `true` iff \f$a\f$ is a primitive element modulo \f$p^e\f$.
         * This method uses the prime factor decomposition of \f$p-1\f$ and its
         * inverse factors.
         */
        bool isPrimitiveElement (const Int & a) const;

        /**
         * Returns `true` iff \f$(-1)^{k+1}V[k]\f$ is a primitive element modulo
         * \f$p^e\f$. This method uses the prime factor decomposition of
         * \f$p-1\f$ and its inverse factors.
         */
        bool isPrimitiveElement (const IntVec &V, int64_t k) const ;

        /**
         * Sets the value of \f$p\f$, \f$e\f$ and \f$m = p^e\f$.
         */
        void setpe (const Int & p, long e);

        /**
         * Gets the value of \f$p\f$.
         */
        Int getP () { return m_p; }

        /**
         * Gets the value of \f$e\f$.
         */
        long getE () { return m_e; }

        /**
         * Sets the value of the factorization `f`.
         */
        void setF (const IntFactorization<Int> & f) { m_f = f; }

        /**
         * Gets the value of the factorization `f`.
         */
        IntFactorization<Int> getF () { return m_f; }

        /**
         * Returns the representation of this object as a string.
         */
        std::string toString () const;

      private:

        /**
         * The prime number \f$p\f$ in this object.
         */
        Int m_p;

        /**
         * Exponent used to compute the modulus \f$m = p^e\f$.
         */
        long m_e;

        /**
         * The modulus \f$m = p^e\f$.
         */
        Int m_m;

        /**
         * The factorization of \f$p-1\f$.
         */
        IntFactorization<Int> m_f;
    };


  template<typename Int>
    PrimitiveInt<Int>::PrimitiveInt () : m_e(1)  {
    m_p = 0;
    m_m = 0;
  }


  //===========================================================================

  template<typename Int>
    void PrimitiveInt<Int>::setpe (const Int & p, long e)  {
      m_p = p;
      m_e = e;
      m_m = NTL::power (m_p, m_e);
    }

  //===========================================================================

  template<typename Int>
    PrimitiveInt<Int>::PrimitiveInt (const IntFactorization<Int> & f,
        const Int & p, long e) : m_f(f)  {
    setpe(p, e);
  }

  //===========================================================================

  template<typename Int>
    std::string PrimitiveInt<Int>::toString () const  {
      std::ostringstream out;
      out << "p = " << m_p << std::endl;
      out << "e = " << m_e << std::endl;
      out << "m = " << m_m << std::endl;
      out << "\nf is the factorization of  " << m_f.toString () << std::endl;
      return out.str ();
    }

  //===========================================================================

  template<typename Int>
    bool PrimitiveInt<Int>::isPrimitiveElement (const Int & a) const  {
      if (0 == m_p)
        throw std::range_error("PrimitiveInt::isPrimitiveElement: p = 0");
      if (0 == a)
        return false;
      Int t1, t2;
      t1 = a;
      if (t1 < 0)
        t1 += m_m;
      const std::vector<Int> invList = m_f.getInvFactorList();
      //   assert (!(invList.empty ()));
      for (auto it = invList.begin(); it != invList.end(); it++) {
        if (*it == (m_m-1)) continue;
        t2 = NTL::PowerMod (t1, *it, m_m);
        if (t2 == 1)
          return false;
      }
      return true;
    }


  //===========================================================================

  template<typename Int>
    bool PrimitiveInt<Int>::isPrimitiveElement (const IntVec & V, int64_t k) const  {
      Int a;
      if (k & 1)
        a = V[k];
      else
        a = -V[k];
      return isPrimitiveElement (a);
    }
  template class PrimitiveInt<std::int64_t>;
  template class PrimitiveInt<NTL::ZZ>;

}
#endif

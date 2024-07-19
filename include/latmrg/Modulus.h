#ifndef LATMRG_MODULUS_H
#define LATMRG_MODULUS_H

#include "latticetester/Util.h"

#include <cassert>


namespace LatMRG {

  /**
   * This class offers a few tools to work with a modulus `m`.
   * It keeps the value of `m` as an `Int` and its square root.
   * It can also verify the max period conditions and compute the reduced modulus
   * for the lattice structure over a single cycle, when `m` is a power of a prime.
   * In the latter case, the modulus must be initialized as \f$m =b^e + c\f$.
   *
   * I am not sure if we should keep this class.  It seems that what it contains should
   * be move to where it is used.    ********
   */
  template<typename Int>
    class Modulus {
      public:

        Modulus ();

        /**
         * Constructor with modulus of congruence \f$m\f$.
         */
        Modulus (const Int & m);

        /**
         * Constructor with value \f$m =b^e + c\f$. Restrictions: \f$b>1\f$ and
         * \f$e > 0\f$.
         */
        Modulus (long b, long e, long c);

        /**
         * Destructor.
         */
        virtual ~Modulus ();

        /**
         * Initializes with value \f$m\f$. Computes `mRac` and `mRacNeg`.
         */
        void init (const Int & m);

        /**
         * Initializes with value \f$m =b^e + c\f$. Restrictions: \f$b>1\f$ and
         * \f$e > 0\f$. Computes `mRac` and `mRacNeg`.
         */
        void init (long b, long e, long c);

        /**
         * Reduces the modulus \f$m\f$ and sets the variable `mRed` to the
         * reduced modulus. The modulus must have the form \f$m=p^e\f$. The
         * multiplier of the LCG is \f$a\f$.
         */
        void reduceM (const Int & a);

        /**
         * Assumes that \f$m\f$ is a power of a prime \f$p=b\f$ and that we have
         * a multiplicative LCG (\f$k = 1\f$).
         * Returns `true` iff the maximal period conditions are satisfied for this `a`.
         * For `b=2`, this holds iff `a mod 8 = 3` or `5`.
         */
        bool perMaxPowPrime (const Int & a);

        /**
         * Value \f$m\f$ of the modulus.
         */
        Int m;

        /**
         * Reduced value of the modulus, in case m is a power of a prime. Computed by `reduceM`.
         */
        Int mRed;

        /**
         * This flag is `true` when \f$m\f$ is prime, otherwise `false`.
         */
        bool primeF;

        /**
         * When this flag is `true`, the value of \f$m\f$ is built out of the
         * three numbers \f$b\f$, \f$e\f$ and \f$c\f$ as described below;
         * otherwise, the flag is set `false`.
         */
        bool threeF;
        long b;
        long e;

        /**
         * When `threeF` is `true`, then \f$m\f$ is given in the form \f$m = b^e +
         * c\f$; otherwise, \f$b\f$, \f$e\f$ and \f$c\f$ are undefined.
         */
        long c;

        /**
         * \f$\sqrt{\lfloor m \rfloor}\f$.
         */
        Int mRac;

        /**
         * \f$-\sqrt{\lfloor m \rfloor}\f$.
         */
        Int mRacNeg;

      private:

        /**
         * The constant \f$b - 1\f$.
         */
        Int bm1;

        /**
         * The constant \f$b^2\f$.
         */
        Int b2;

        /**
         * Work variables.
         */
        Int Y;  // b as an Int.
        Int Eight(8);
        Int Four(4);

    }; // End class declaration

  //===========================================================================

  template<typename Int>
    Modulus<Int>::Modulus()
    {}


  //===========================================================================

  template<typename Int>
    Modulus<Int>::Modulus(const Int & j)
    {
      init(j);
    }


  //===========================================================================

  template<typename Int>
    Modulus<Int>::Modulus (long b, long e, long c)   {
      init(b, e, c);
    }


  //===========================================================================

  template<typename Int>
    Modulus<Int>::~Modulus()
    {}


  //===========================================================================

  template<typename Int>
    void Modulus<Int>::init (const Int & j)
    {
      mRed = m = j;
      threeF = false;
      b = 0;
      e = 0;
      c = 0;
      mRac = (Int) NTL::SqrRoot (m);
      mRacNeg = -mRac;
    }


  //===========================================================================

  template<typename Int>
    void Modulus<Int>::init (long b, long e, long c)   {
      assert (b > 1);
      m = b;
      if (e <= 0) {
        init (m);
        return;
      }
      m = NTL::power (m, e) + c;
      init (m);
      this->b = b;
      this->e = e;
      this->c = c;
      threeF = true;
      bm1 = b - 1;
      NTL::conv (Y, b);
      b2 = Y * Y;
    }


  //===========================================================================

  template<typename Int>
    bool Modulus<Int>::perMaxPowPrime (const Int & a)   {
      if (!threeF || c != 0) {
        LatticeTester::MyExit(1, "perMaxPowPrime:   m must be a power of a prime");
        return false;
      }
      if (b == 2) {
        // LatticeTester::Modulo (a, Eight, Y);
        Y = a % Eight;
        return (Y == 5) || (Y == 3);
      } else {
        Y = a % b2;
        Y = NTL::PowerMod (Y, bm1, b2);
        return (Y != 1);
      }
    }

  //===========================================================================

  template<typename Int>
    void Modulus<Int>::reduceM (const Int & a)
    {
      if (!threeF || c != 0)
        return;
      Int Y2, Y3, Y4;

      // We now assume that m is a power of a prime b
      if (b == 2) {
        // m is a power of 2
        LatticeTester::Modulo (a, Eight, Y);
        if (Y == 5)
          LatticeTester::Quotient (m, Four, mRed);
        else if (Y == 3)
          LatticeTester::Quotient (m, Eight, mRed);
        else {
          Y2 = a;
          Y3 = Eight;
          if (Y == 1)
            --Y2;
          else if (Y == 7)
            ++Y2;
          else
            LatticeTester::MyExit (1, "Power of 2 modulus with even multiplier");
          do {
            Y3 = Y3 + Y3;
            LatticeTester::Modulo (Y2, Y3, Y4);
          } while (Y4 == 0);
          LatticeTester::Quotient (m, Y3, mRed);
        }

      } else {
        // mj is a power of a prime b > 2
        NTL::conv (Y3, b);
        if (perMaxPowPrime (a))
          LatticeTester::Divide (mRed, Y4, m, Y3);
        else {
          Y2 = b2;
          do {
            Y2 *= b;        // Y = b^e
            Y = a % Y2;
            Y = NTL::PowerMod (Y, bm1, Y2);
          } while (Y == 1);
          LatticeTester::Divide (Y2, Y4, Y2, Y3);
          LatticeTester::Quotient (m, Y2, mRed);
        }
      }
    }

} // End namespace LatMRG
#endif

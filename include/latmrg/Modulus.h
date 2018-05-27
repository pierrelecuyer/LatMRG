#ifndef MODULUS_H
#define MODULUS_H
#include "latticetester/Util.h"

#include <cassert>


namespace LatMRG {

  /**
   * This class keeps parameters closely associated with a modulus of
   * congruence. Using it, it will not be necessary to recalculate the square
   * roots of large integers, which are used repeatedly in searches for good
   * generators.
   *
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
         * Assumes that \f$m\f$ is a power of a prime \f$p=b\f$, the order \f$k
         * = 1\f$, and the recurrence is homogeneous. Returns `true` iff the
         * maximal period conditions are satisfied.
         */
        bool perMaxPowPrime (const Int & a);

        /**
         * Value \f$m\f$ of the modulus.
         */
        Int m;

        /**
         * Reduced value of the modulus. Computed by `reduceM`.
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
        Int Y, Eight, Four;
    }; // End class declaration

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
    Modulus<Int>::Modulus (long m1, long m2, long m3)
    {
      init(m1, m2, m3);
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
      mRac = (Int) SqrRoot (m);
      mRacNeg = -mRac;
      Eight = 8;
      Four = 4;
    }


  //===========================================================================

  template<typename Int>
    void Modulus<Int>::init (long m1, long m2, long m3)
    {
      assert (m1 > 1);
      m = m1;
      if (m2 <= 0) {
        init (m);
        return;
      }

      m = power (m, m2) + m3;
      init (m);
      b = m1;
      e = m2;
      c = m3;
      threeF = true;
      bm1 = b - 1;
      conv (Y, b);
      b2 = Y*Y;
    }


  //===========================================================================

  template<typename Int>
    bool Modulus<Int>::perMaxPowPrime (const Int & A)
    {
      if (!threeF || c != 0) {
        LatticeTester::MyExit(1, "perMaxPowPrime:   m must be a power of a prime");
        return false;
      }

      if (b == 2) {
        LatticeTester::Modulo (A, Eight, Y);
        return (Y == 5) || (Y == 3);
      } else {
        Y = A % b2;
        Y = PowerMod (Y, bm1, b2);
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
        conv (Y3, b);
        if (perMaxPowPrime (a))
          LatticeTester::Divide (mRed, Y4, m, Y3);
        else {
          Y2 = b2;
          do {
            Y2 *= b;        // Y = b^e
            Y = a % Y2;
            Y = PowerMod (Y, bm1, Y2);
          } while (Y == 1);
          LatticeTester::Divide (Y2, Y4, Y2, Y3);
          LatticeTester::Quotient (m, Y2, mRed);
        }
      }
    }



} // End namespace LatMRG
#endif

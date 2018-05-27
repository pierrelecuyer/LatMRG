#ifndef MRGCOMPONENTFACTORY_H
#define MRGCOMPONENTFACTORY_H
#include "latmrg/MRGComponent.h"
#include "latticetester/Util.h"


namespace LatMRG {

  /**
   * This class is used to create <tt>MRGComponent</tt>s from other types of
   * recurrences [for instance, multiply-with-carry (MWC)].
   *
   */
  template<typename Int>
    class MRGComponentFactory {
      public:

        /**
         * Creates a `MRGComponent` from an multiply-with-carry type of 
         * recurrence. The MWC recurrence has the form
         * \f{align*}{
         *    x_n
         *    &
         *    =
         *    (a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1})d \mbox{ mod } b,
         *  \\
         *   c_n
         *    &
         *    =
         *  \lfloor(a_0 x_n + a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1}) / b \rfloor
         * \f}
         * where \f$b\f$ is a positive integer, \f$a_0,…,a_k\f$ are arbitrary
         * integers such that \f$a_0\f$ is relatively prime to \f$b\f$, and 
         * \f$d\f$ is the multiplicative inverse of \f$-a_0\f$ modulo \f$b\f$. 
         * The MRG derived from such a MWC is defined by
         * \f[
         *   m = \sum^k_{i=0} a_i b^i
         * \f]
         * where \f$a\f$ is the inverse of \f$b\f$ in arithmetic modulo 
         * \f$m\f$.
         */
        static MRGComponent<Int> * fromMWC (const Int & b, const MVect & a, int k);
    }; // End class declaration

  template<typename Int>
    MRGComponent<Int>* MRGComponentFactory<Int>::fromMWC(const Int & b,
        const MVect & a, int k)
    {
      Int _m;
      Int _b;
      MVect _a;

      Int d, e, f, g;

      conv(_m, 0);
      conv(_b, 1);
      LatticeTester::CreateVect(_a, 1);

      //Calcul de m
      for (int i = 0; i <= k; i++) {
        _m += a[i] * _b;
        _b *= b;
      }

      //Calcul de a
      LatticeTester::Euclide(b, _m, _a[1], d, e, f, g);

      if (g < 0)
        _a[1] = -_a[1];

      MRGComponent<Int>* lat = new MRGComponent<Int>(_m, _a, 1);
      LatticeTester::DeleteVect(_a);
      return lat;
    }

} // End namespace LatMRG
#endif

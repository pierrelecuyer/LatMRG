#ifndef MRGLATTICEFACTORY_H
#define MRGLATTICEFACTORY_H
#include "latticetester/Const.h"
#include "latmrg/Const.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"


namespace LatMRG {

  /**
   * This class is used to create <tt>MRGLattice</tt>’s from other types of
   * recurrences than a MRG or from a combination of several MRG lattices.
   *
   */
  template<typename Int, typename Dbl>
    class MRGLatticeFactory {
      private:
        typedef NTL::vector<Int> IntVec;
      public:

        /**
         * Creates a `MRGLattice` from the combination of \f$J\f$
         * <tt>MRGComponent</tt>’s, with maximal dimension `maxDim`, lacunary 
         * indices `Lac`, lattice type `lat` and vector norm `norm`. `Lac` can 
         * be `NULL` if no lacunary indices are to be used. The combined MRG is 
         * calculated as described in \cite rLEC96b.
         */
        static MRGLattice<Int, Dbl> * fromCombMRG (MRGComponent<Int> **comp, int J,
            int maxDim, IntVec * Lac, LatticeType lat,
            LatticeTester::NormType norm);

        /**
         * Creates a `MRGLattice` from a multiply-with-carry (MWC) recurrence,
         * with maximal dimension `maxDim`, lacunary indices `Lac`, lattice
         * type `lat` and vector norm `norm`. The MWC recurrence is
         * \f{align*}{
         *    x_n
         *    &
         *    =
         *    (a_1 x_{n-1} + \cdots+ a_k x_{n-k} + c_{n-1})d \mbox{ mod } b,
         *  \\
         *   c_n
         *    &
         *    =
         *  \lfloor(a_0 x_n + a_1 x_{n-1} +\cdots+ a_k x_{n-k} + c_{n-1}) / b \rfloor
         * \f}
         * where \f$b\f$ is a positive integer, \f$a_0,\cdots,a_k\f$ are 
         * arbitrary integers such that \f$a_0\f$ is relatively prime to 
         * \f$b\f$, and \f$d\f$ is the multiplicative inverse of \f$-a_0\f$ mod 
         * \f$b\f$. The MRG derived from such a MWC is defined by 
         * \f$ m = \sum^k_{l=0} a_l b^l \f$ and \f$a\f$ is the inverse of 
         * \f$b\f$ in arithmetic modulo \f$m\f$.
         */
        static MRGLattice<Int, Dbl> * fromMWC (const IntVec & a, const Int & b, 
            int maxDim, int k, IntVec *Lac, LatticeType lat,
            LatticeTester::NormType norm);

        /**
         * Same as above, but with no lacunary indices.
         */
        static MRGLattice<Int, Dbl> * fromMWC (const IntVec & a, const Int & b,
            int maxDim, int k, LatticeType lat,
            LatticeTester::NormType norm);
    }; // End class declaration

  //===========================================================================

  template<typename Int, typename Dbl>
    MRGLattice<Int, Dbl> *MRGLatticeFactory<Int, Dbl>::fromCombMRG (MRGComponent<Int> ** comp,
        int J, int maxDim, IntVec * I, LatticeType type,
        LatticeTester::NormType norm)
    {
      Int _m;
      int _k = 1;

      /*
Remark: before the code was
''
Int _n[J];

''
But this sends a compiler error: "variable length array of non-POD
type Int".
So instead we decied to use a vector to perform the same operations
with _n[].
*/
      IntVec _n;
      _n.resize(J);

      Int d, e, f, g;

      NTL::conv (_m, 1);
      for (int j = 0; j < J; j++) {
        _k = std::max (_k, comp[j]->k);
        _m *= comp[j]->getM();
      }

      // Calcul de nj
      Int tmp;
      for (int j = 0; j < J; j++) {
        LatticeTester::Quotient (_m, comp[j]->getM(), tmp);
        LatticeTester::Euclide (tmp, comp[j]->getM(), _n[j], d, e, f, g);
        if (g < 0)
          _n[j] = -_n[j];
        _n[j] *= tmp;
        comp[j]->nj = _n[j];
      }

      IntVec _a;
      LatticeTester::CreateVect (_a, _k);

      // Calcul de ai
      for (int i = 1; i <= _k; ++i) {
        NTL::clear (_a[i]);
        for (int j = 0; j < J; ++j) {
          if (comp[j]->k >= i) {
            tmp = comp[j]->a[i] * _n[j];
            _a[i] += tmp;
          }
        }
        LatticeTester::Modulo (_a[i], _m, _a[i]);
      }
      for (int j = 0; j < J; ++j) {
        tmp = _a[1] - comp[j]->a[1];
        LatticeTester::Modulo (tmp, comp[j]->getM(), tmp);
        assert (0 == tmp);
      }

      // Calcul de rho, lossrho et rhoj
      Int rho;
      Int lossRho;
      rho = 1;
      lossRho = 1;

      for (int j = 0; j < J; j++) {
        comp[j]->rho = NTL::power (comp[j]->getM(), comp[j]->k) - 1;
        //      LatticeTester::Euclide (g, rho, tmp, d, e, f, comp[j]->rho);
        LatticeTester::Euclide (comp[j]->rho, rho, tmp, d, e, f, g);
        lossRho *= g;
        rho *= comp[j]->rho;
        LatticeTester::Quotient (rho, g, rho);
      }

      MRGLattice<Int, Dbl> *lat;

      if (I != 0)
        lat = new MRGLattice<Int, Dbl> (_m, _a, maxDim, _k, *I, type, norm);
      else
        lat = new MRGLattice<Int, Dbl> (_m, _a, maxDim, _k, type, norm);

      lat->setRho (rho);
      lat->setLossRho (lossRho);
      LatticeTester::DeleteVect (_a);
      lat->comp.clear();
      lat->comp.reserve(J);
      MRGComponent<Int> *mycomp;
      for (int j = 0; j < J; j++) {
        mycomp = new MRGComponent<Int> (*comp[j]);
        lat->comp.push_back(mycomp);
      }
      return lat;
    }


  //===========================================================================
#if 0

  // maybe an old and not updated implementation ?

  MRGLattice *MRGLatticeFactory::fromCombMRG (const Int * m,
      const MMat & coef, int MaxDim, int *k, int J, BVect * I, LatticeType lat,
      NormType norm)
  {
    Int _m;
    int _k = 0;
    Int _n[J];
    IntVec _a;

    Int d, e, f, g;

    NTL::conv (_m, 1);

    // Calcul de m et k;
    for (int i = 0; i < J; ++i) {
      _k = std::max (_k, k[i]);
      _m *= m[i];
    }

    // Calcul de nj
    Int t1;
    for (int i = 0; i < J; ++i) {
      t1 = _m / m[i];

      Euclide (t1, m[i], _n[i], d, e, f, g);

      if (g < 0) {
        _n[i] = -_n[i];
      }
      _n[i] *= t1;

    }

    CreateVect (_a, _k);
    Int tmp;
    // Calcul de ai
    for (int i = 1; i <= _k; ++i) {
      clear (_a[i]);
      for (int j = 0; j < J; ++j) {
        if (k[j] >= i) {
          tmp = coef[j][i] * _n[j];
          _a[i] += tmp;
        }
      }
      Modulo (_a[i], _m, _a[i]);
    }

    if (I != 0) {
      return new MRGLattice (_m, _a, MaxDim, _k, *I, lat, norm);
    } else {
      return new MRGLattice (_m, _a, MaxDim, _k, lat, norm);
    }
  }


  //=========================================================================

  MRGLattice *MRGLatticeFactory::fromCombMRG (const Int m[],
      const MMat & coef, int MaxDim, int k[], int J, LatticeType lat,
      NormType norm)
  {
    return fromCombMRG (m, coef, MaxDim, k, J, 0, lat, norm);
  }
#endif


  //===========================================================================

  template<typename Int, typename Dbl>
    MRGLattice<Int, Dbl> *MRGLatticeFactory<Int, Dbl>::fromMWC (const IntVec & a,
        const Int & b, int MaxDim, int k, LatticeType lat,
        LatticeTester::NormType norm)
    {
      return fromMWC (a, b, MaxDim, k, 0, lat, norm);
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    MRGLattice<Int, Dbl> *MRGLatticeFactory<Int, Dbl>::fromMWC (const IntVec & a,
        const Int & b, int MaxDim, int k, IntVec * I, LatticeType lat_t,
        LatticeTester::NormType norm)
    {
      Int _m;
      Int _b;
      IntVec _a;

      Int d, e, f, g;

      NTL::conv (_m, 0);
      NTL::conv (_b, 1);
      LatticeTester::CreateVect (_a, 1);

      // Calcul de m
      for (int i = 0; i <= k; i++) {
        _m += a[i] * _b;
        _b *= b;
      }

      // Calcul de a
      LatticeTester::Euclide (b, _m, _a[1], d, e, f, g);

      if (g < 0) {
        _a[1] = -_a[1];
      }

      MRGLattice<Int, Dbl> *lat;
      if (I) {
        lat = new MRGLattice<Int, Dbl> (_m, _a, MaxDim, 1, *I, lat_t, norm);
      } else {
        lat = new MRGLattice<Int, Dbl> (_m, _a, MaxDim, 1, lat_t, norm);
      }
      LatticeTester::DeleteVect (_a);
      return lat;
    }

  extern template class MRGLatticeFactory<std::int64_t, double>;
  extern template class MRGLatticeFactory<NTL::ZZ, double>;
  extern template class MRGLatticeFactory<NTL::ZZ, NTL::RR>;

} // End namespace LatMRG
#endif

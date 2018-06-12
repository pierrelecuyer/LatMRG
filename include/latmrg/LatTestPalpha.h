#ifndef LATTESTPALPHA_H
#define LATTESTPALPHA_H
#include "latticetester/Normalizer.h"

#include "latmrg/LatConfig.h"
#include "latmrg/LatticeTest.h"


namespace LatMRG {

  /**
   * This class applies the \f$P_{\alpha}\f$ test by implementing the abstract
   * class `LatticeTest`. The figure of merit is the
   * \f$P_{\alpha}\f$ criterion (see class `Palpha` on page (FIXME: page#)).
   * The main program is obtained by compiling the `LatMain.cc` file, and the
   * executable file is called `latLLDD`, since the test is implemented only
   * for a small number of points (\f$m < 2^{31}\f$). See the description of
   * the `LatMain` program on page (FIXME: page#).
   *
   */
  template<typename Int>
    class LatTestPalpha : public LatticeTest {
      public:

        /**
         * Constructor. The \f$P_{\alpha}\f$ test will be applied on lattice whose
         * parameters are in `config`. The `bounds` \f$B_{\alpha}(s)\f$ may be used
         * to normalize the \f$P_{\alpha}(s)\f$ values.
         */
        LatTestPalpha (LatticeTester::Normalizer<RScal> * bounds, LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal> * lat);

        /**
         * Destructor.
         */
        ~LatTestPalpha () {};
        void setConfig (LatConfig<Int> * config) {m_config = config;};

        /**
         * Applies the \f$P_{\alpha}\f$ test for dimensions varying from
         * <tt>fromDim</tt> to `toDim`. Whenever the normalized value of the
         * merit is smaller than `minVal` for any dimension, the method returns
         * `false` immediately. The method returns `false` if the test was
         * interrupted for any reason before completion, and it returns `true`
         * upon success. The results of the last test are kept in `merit`.
         */
        bool test (int fromDim, int toDim, double minVal[]);

      private:

        /**
         * Contains the parameters of the test.
         */
        LatConfig<Int> *m_config;

        /**
         * The \f$B_{\alpha}\f$ bounds used to normalize the
         * \f$P_{\alpha}\f$.
         */
        LatticeTester::Normalizer<RScal> *m_bound;

        /**
         * Prepares and dispatches the results for dimension `dim` to all
         * observers attached to this test.
         */
        void prepAndDisp (int dim);

        /**
         * Sends the initialization message to all observers attached to this
         * test.
         */
        void init ();
    }; // End class declaration

  //===========================================================================

  template<typename Int>
    void LatTestPalpha<Int>::init ()
    {
      m_merit.setDim(m_toDim);
      const int N = 3;
      std::string header[N];
      if (m_config->calcPalpha == LatticeTester::NORMPAL) {
        header[0] = "P_a";
        header[1] = "P_a/B_a";
        header[2] = "Cumul CPU t(sec)";
        dispatchTestInit ("Palpha", header, N);
      } else {
        header[0] = "P_a";
        header[1] = "Cumul CPU t(sec)";
        dispatchTestInit ("Palpha", header, N - 1);
      }
      timer.init();
    }


  //===========================================================================

  template<typename Int>
    void LatTestPalpha<Int>::prepAndDisp (int dim)
    {
      const int N = 3;
      double results[N];
      if (m_config->calcPalpha == LatticeTester::NORMPAL) {
        results[0] = m_merit.getMerit (dim);
        results[1] = m_merit[dim];
        results[2] = timer.val(Chrono::SEC);
        dispatchResultUpdate (results, N);
      } else {
        results[0] = m_merit[dim];
        results[1] = timer.val(Chrono::SEC);
        dispatchResultUpdate (results, N - 1);
      }
    }


  //===========================================================================

  template<typename Int>
    LatTestPalpha<Int>::LatTestPalpha (LatticeTester::Normalizer<RScal> * normal, LatticeTester::IntLattice<Int, BScal, BVect, BMat, NScal, NVect, RScal> * lat): LatticeTest (lat)
  {
    m_criter = LatticeTester::PALPHA;
    m_bound = normal;
    //   m_config = config;
    //   m_fromDim = m_config->td[0];
    //   m_toDim = m_config->td[1];
    m_dualF = false;
    m_maxAllDimFlag = false;
  }


  //===========================================================================
  /*
     bool LatTestPalpha::test (double MinVal[])
     {
     return test (m_config->fromDim, m_config->toDim, MinVal);
     }
     */

  //===========================================================================

  template<typename Int>
    bool LatTestPalpha<Int>::test (int fromDim, int toDim, double minVal[])
    {
      m_fromDim = fromDim;
      m_toDim = toDim;
      init ();

      if (m_config->verifyM) {
        if (LatticeTester::PRIME == LatticeTester::IntFactor<Int>::isPrime (m_config->comp[0]->getM(), 50)) {
          std::cout << "Verify prime m:   true" << std::endl;
          m_config->primeM = true;
        } else {
          std::cout << "Verify prime m:   false" << std::endl;
          m_config->primeM = false;
        }
      }

      if (m_config->verifyP) {
        MRGComponent<Int> mrg (m_config->comp[0]->getM(), 1, DECOMP, 0, DECOMP_PRIME,  0);
        if (mrg.maxPeriod (m_config->comp[0]->a)) {
          std::cout << "Verify maximal period:   true" << std::endl;
          m_config->maxPeriod = true;
        } else {
          std::cout << "Verify maximal period:   false" << std::endl;
          m_config->maxPeriod = false;
        }
        std::cout << std::endl;
      }

      //PW_TODO pourquoi c'est commenté ?
      //PalphaLCG palpha (*m_config);
      //const int alpha = m_config->alpha;
      //double x;
      //int dim;
      /*
         if (m_config->calcPalpha == PAL || m_config->calcPalpha == NORMPAL) {

      // SI LE MODULO EST PREMIER  // Si periode maximale
      if (m_config->primeM && m_config->maxPeriod) {
      for (dim = fromDim; dim <= toDim; dim++) {
      switch (alpha) {
      case 2:
      x = palpha.calcPalpha2 (dim);
      break;
      case 4:
      x = palpha.calcPalpha4 (dim);
      break;
      case 6:
      x = palpha.calcPalpha6 (dim);
      break;
      case 8:
      x = palpha.calcPalpha8 (dim);
      break;
      default:
      cerr << " Valeur de alpha invalide" << std::endl;
      return false;
      }

      if (m_config->calcPalpha == PAL)
      m_merit[dim] = x;
      else {
      m_merit.getMerit(dim) = x;
      if (m_bound->getCst(dim) < 0.0)
      m_merit[dim] = -1.0;
      else
      m_merit[dim] = x / m_bound->getCst(dim);
      }
      prepAndDisp (dim);
      }

      } else {
      for (dim = fromDim; dim <= toDim; dim++) {
      switch (alpha) {
      case 2:
      x = palpha.calcPalpha2PerNonMax (dim);
      break;
      case 4:
      x = palpha.calcPalpha4PerNonMax (dim);
      break;
      case 6:
      x = palpha.calcPalpha6PerNonMax (dim);
      break;
      case 8:
      x = palpha.calcPalpha8PerNonMax (dim);
      break;
      default:
      cerr << " Valeur de alpha invalide" << std::endl;
      return false;
      }
      if (m_config->calcPalpha == PAL)
      m_merit[dim] = x;
      else {
      m_merit.getMerit(dim) = x;
      if (m_bound->getCst(dim) < 0.0)
      m_merit[dim] = -1.0;
      else
      m_merit[dim] = x / m_bound->getCst(dim);
      }
      prepAndDisp (dim);
      }
      }

      } else if (m_config->calcPalpha == BAL) {
      for (dim = fromDim; dim <= toDim; dim++) {
      x = m_bound->getCst(dim);
      m_merit[dim] = x;
      prepAndDisp (dim);
    }

    }
    */
      return true;
    }

} // End namespace LatMRG
#endif

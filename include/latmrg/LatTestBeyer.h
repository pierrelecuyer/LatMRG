#ifndef LATTESTBEYER_H
#define LATTESTBEYER_H

#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"

#include "latmrg/LatticeTest.h"


namespace LatMRG {

  /**
   * This class implements the *Beyer* test. It implements the abstract class
   * `LatticeTest`. The figure of merit for this test is the Beyer quotient.
   * The main program is obtained by compiling the `LatMain.cc` file (see the
   * description of module `LatMain`, on page (FIXME: page#), for how to run a
   * program).
   * \remark **Richard:** Dixit Pierre: il est douteux que cette classe devrait
   * exister: ce devrait être une méthode de Lattice ou quelque chose du genre.
   * Idem pour les autres `LatTest*`
   *
   */
  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
      class LatTestBeyer : public LatticeTest {
        public:

          /**
           * Constructor. The *Beyer* test will be applied on lattice `lat`.
           */
          LatTestBeyer (LatticeTester::IntLattice<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl> * lat);

          /**
           * Destructor.
           */
          ~LatTestBeyer () {};

          /**
           * Applies the *Beyer* test for dimensions varying from
           * <tt>fromDim</tt> to `toDim`. Whenever the normalized value of the
           * merit is smaller than `minVal` for any dimension, the method returns
           * `false` immediately. The method returns `false` if the test was
           * interrupted for any reason before completion, and it returns `true`
           * upon success. The results of the last test are kept in `merit`.
           */
          bool test (int fromDim, int toDim, double minVal[]);

        private:

          /**
           * Prepares and dispatches the results for dimension `dim` to all observers
           * attached to this test.
           */
          void prepAndDisp (int dim);

          /**
           * Sends the initialization message to all observers attached to this
           * test.
           */
          void init ();
      }; // End class LatTestBeyer

  // Implementation of class LatTestBeyer
  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
      LatTestBeyer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl>::LatTestBeyer (LatticeTester::IntLattice<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl> * lat): LatticeTest (lat)
    {
      m_criter = LatticeTester::BEYER;
    }

  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
      bool LatTestBeyer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl>::test (int fromDim, int toDim, double minVal[])
      {
        m_fromDim = fromDim;
        m_toDim = toDim;
        init ();

        resetFromDim (m_lat->getOrder (), fromDim);
        while (m_lat->getDim () < fromDim)
          m_lat->incDim ();
        LatticeTester::Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RVect, RMat> red (*m_lat);

        m_lat->dualize ();
        red.preRedDieter (0);
        m_lat->dualize ();

        while (true) {
          if (m_dualF)
            m_lat->dualize ();

          //cout << "CURRENT BASIS =\n";
          //testPrinting(m_lat->getBasis(), "testMatrix");

          bool success = red.reductMinkowski (0);
          int dim = m_lat->getDim ();
          if (success) {
            m_lat->updateScalL2Norm (0);
            m_lat->updateScalL2Norm (dim-1);

            double x1, x2;        // si VV[1] et VV[dim] sont tres
            // grands, il faudrait envisager de changer x1 et x2 en xdouble.
            NTL::conv (x1, m_lat->getVecNorm (0));
            NTL::conv (x2, m_lat->getVecNorm (dim-1));
            m_merit[dim] = x1 / x2;


            /*
               cout << "\nlac *** dim = " << dim << endl;
               cout << "lac *** x1 = " << x1 << endl;
               cout << "lac *** x2 = " << x2 << endl;
               cout << "lac *** m_merit = " << sqrt(m_merit[dim]) << endl;
               */


            // Si on sait deja que ce gen. ne pourra etre retenu,
            // on le rejette tout de suite et on arrete le test.
            if ((m_maxAllDimFlag && (m_merit[dim] < minVal[toDim]))
                || (m_merit[dim] < minVal[dim])) {
              m_merit[dim] = 0.0;
              return false;
            }
            if (3 == m_detailF) {
              this->dispatchLatUpdate(*m_lat);
            }

            prepAndDisp (dim);
            if (m_dualF)
              m_lat->dualize ();

          } else {
            m_merit[dim] = 0.0;
            return false;
          }

          if (dim == toDim)
            break;
          m_lat->incDim();
          red = LatticeTester::Reducer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl, RVect, RMat>(*m_lat);
        }

        return true;
      }

  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
      void LatTestBeyer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl>::init ()
      {
        m_merit.setDim(m_toDim);
        const int N = 2;
        std::string header[N];
        header[0] = "q_t";
        header[1] = "Cumul CPU t(sec)";
        this->dispatchTestInit ("Beyer", header, N);
        this->timer.init();
      }

  //===========================================================================

  template<typename Int, typename BasInt, typename BasIntVec,
    typename BasIntMat, typename Dbl, typename DblVec, typename RedDbl>
      void LatTestBeyer<Int, BasInt, BasIntVec, BasIntMat, Dbl, DblVec, RedDbl>::prepAndDisp (int dim)
      {
        const int N = 2;
        double results[N];
        results[0] = sqrt (m_merit[dim]);
        results[1] = this->timer.val (Chrono::SEC);
        dispatchResultUpdate (results, N);
      }

} // End namespace LatMRG

#endif // LATTESTBEYER_H

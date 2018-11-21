#ifndef TESTPROJECTIONS_H
#define TESTPROJECTIONS_H

#include <cstring>
#include <sstream>

#include "latticetester/Writer.h"
#include "latticetester/WriterRes.h"
#include "latticetester/UniformWeights.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Types.h"
#include "latticetester/Weights.h"

#include "latmrg/ProjIterator.h"
#include "latmrg/MeritProj.h"
#include "latmrg/LatticeTest.h"
#include "latmrg/ProjIteratorSuccCoords.h"
#include "latmrg/ProjIteratorNonSuccCoords.h"

namespace
{
  const char* espaceM = " \t\t   ";
  const char* espaceG = "                                                   ";

  const LatticeTester::UniformWeights unitWeights(1.0);

  const std::string formatIndices (const LatticeTester::Coordinates & ens)
  {
    std::ostringstream os;
    os << "  ";
    LatticeTester::Coordinates::const_iterator it = ens.begin ();
    os << *it;
    ++it;
    while (it != ens.end ()) {
      os << "," << *it;
      ++it;
    }
    os << espaceM;
    return os.str ();
  }

  template<typename Int>
    void printMerit (bool racF, double merit, LatticeTester::Writer<Int> * rw)
    {
      if (racF) {
        rw->writeString (espaceG);
        rw->writeString ("min merit = ");
        rw->writeDouble (LatticeTester::mysqrt (merit));
      } else {
        rw->writeString (espaceG);
        rw->writeString ("min merit = ");
        rw->writeDouble (merit);
      }
      rw->newLine ();
      rw->newLine ();
    }

  // Prints `len` and `merit` on the output of `rw`.
  template<typename Int>
    void printMerLen (bool invF, bool racF, double len, double merit,
        LatticeTester::Writer<Int> * rw)
    {
      double x = len;
      if (racF)
        x = LatticeTester::mysqrt (len);
      if (invF)
        x = 1.0 / x;

      if (racF) {
        rw->writeDouble (x);
        rw->writeString (espaceM);
        rw->writeDouble (LatticeTester::mysqrt (merit));
      } else {
        rw->writeDouble (x);
        rw->writeString (espaceM);
        rw->writeDouble (merit);
      }
      rw->newLine ();
    }

}                                 // namespace

namespace LatMRG {

  /**
   * Implements methods used to calculate the worst-case figure of merit,
   * defined as follows:
   * \f[
   *   M_{t_1,…,t_d} = \min\left[ \min_{k+1\le t\le t_1} \frac{\ell_t}{\ell_t^*(m^k)},\; \min_{2\le s\le k}\; \min_{I\in S(s,t_s)} \frac{\ell_I}{m}, \; \min_{k+1\le s\le d} \;\min_{I\in S(s,t_s)} \frac{\ell_I}{\ell_s^*(m^k)} \right],
   * \f]
   * where \f$S(s, t_s) = \left\{I = \{i_1,…,i_s\} \mid1 = i_1 < \cdots< i_s
   * \le t_s\right\}\f$. In other words, this figure of merit applies the
   * chosen test over the \f$s\f$ successive dimensions for all \f$s
   * \le t_1\f$, and over the \f$r\f$ nonsuccessive dimensions of the lattice
   * of dimension \f$t_r\f$ for all \f$2 \le r \le d\f$. For *dimension
   * stationary* lattices, for example Korobov lattices, only the sets of
   * dimensions whose first coordinate is 1 need to be considered because being
   * dimension stationary means that shifting the lattice by a constant does not
   * change the projection.
   *
   * Here is an example of the spectral test with different projections for a
   * Korobov lattice with \f$m=1021\f$ and \f$a=333\f$, when we consider the
   * dual lattice with normalization BESTLAT. The figure of merit is
   * \f$M_{10,8,6,5}\f$ and the resulting minimal merit is 0.440728, for
   * projection 1,6. The length is the length of the shortest vector of the
   * lattice with the \f${\mathcal{L}}_2\f$ norm.
   *
   * <center>
   *
   * <table class="LatSoft-table LatSoft-has-hlines">
   * <tr class="bt">
   *   <td class="l">projections</td>
   *   <td class="l">length</td>
   *   <td class="l">merit</td>
   * </tr><tr class="bt">
   *   <td class="l">2</td>
   *   <td class="l">22.2036</td>
   *   <td class="l">0.64666</td>
   * </tr><tr>
   *   <td class="l">3</td>
   *   <td class="l">7</td>
   *   <td class="l">0.619324</td>
   * </tr><tr>
   *   <td class="l">4</td>
   *   <td class="l">3.87298</td>
   *   <td class="l">0.576145</td>
   * </tr><tr>
   *   <td class="l">5</td>
   *   <td class="l">3.74166</td>
   *   <td class="l">0.760239</td>
   * </tr><tr>
   *   <td class="l">6</td>
   *   <td class="l">2</td>
   *   <td class="l">0.488395</td>
   * </tr><tr>
   *   <td class="l">7</td>
   *   <td class="l">2</td>
   *   <td class="l">0.552276</td>
   * </tr><tr>
   *   <td class="l">8</td>
   *   <td class="l">2</td>
   *   <td class="l">0.594822</td>
   * </tr><tr>
   *   <td class="l">9</td>
   *   <td class="l">2</td>
   *   <td class="l">0.654906</td>
   * </tr><tr>
   *   <td class="l">10</td>
   *   <td class="l">2</td>
   *   <td class="l">0.697213</td>
   * </tr><tr class="bt">
   *   <td class="l">1,3</td>
   *   <td class="l">25.4951</td>
   *   <td class="l">0.742522</td>
   * </tr><tr>
  *   <td class="l">1,4</td>
    *   <td class="l">20.6155</td>
    *   <td class="l">0.600409</td>
    * </tr><tr>
    *   <td class="l">1,5</td>
    *   <td class="l">30.6105</td>
    *   <td class="l">0.891502</td>
    * </tr><tr>
    *   <td class="l">1,6</td>
    *   <td class="l">15.1327</td>
    *   <td class="l">0.440728</td>
    * </tr><tr>
    *   <td class="l">1,7</td>
    *   <td class="l">16.6433</td>
    *   <td class="l">0.484722</td>
    * </tr><tr>
    *   <td class="l">1,8</td>
    *   <td class="l">32.3883</td>
    *   <td class="l">0.943279</td>
    * </tr><tr class="bt">
    *   <td class="l">1,2,4</td>
    *   <td class="l">7.07107</td>
    *   <td class="l">0.625612</td>
    * </tr><tr>
    *   <td class="l">1,2,5</td>
    *   <td class="l">8.12404</td>
    *   <td class="l">0.718773</td>
    * </tr><tr>
    *   <td class="l">1,2,6</td>
    *   <td class="l">9.48683</td>
    *   <td class="l">0.839346</td>
    * </tr><tr>
    *   <td class="l">1,3,4</td>
    *   <td class="l">9.43398</td>
    *   <td class="l">0.83467</td>
    * </tr><tr>
    *   <td class="l">1,3,5</td>
    *   <td class="l">9.69536</td>
    *   <td class="l">0.857795</td>
    * </tr><tr>
    *   <td class="l">1,3,6</td>
    *   <td class="l">8.12404</td>
    *   <td class="l">0.718773</td>
    * </tr><tr>
    *   <td class="l">1,4,5</td>
    *   <td class="l">7.54983</td>
    *   <td class="l">0.66797</td>
    * </tr><tr>
    *   <td class="l">1,4,6</td>
    *   <td class="l">8.12404</td>
    *   <td class="l">0.718773</td>
    * </tr><tr>
    *   <td class="l">1,5,6</td>
    *   <td class="l">6.78233</td>
    *   <td class="l">0.600066</td>
    * </tr><tr class="bt">
    *   <td class="l">1,2,3,5</td>
    *   <td class="l">5.74456</td>
    *   <td class="l">0.854561</td>
    * </tr><tr>
    *   <td class="l">1,2,4,5</td>
    *   <td class="l">3.87298</td>
    *   <td class="l">0.576145</td>
    * </tr><tr>
    *   <td class="l">1,3,4,5</td>
    *   <td class="l">5.47723</td>
    *   <td class="l">0.814792</td>
    * </tr>
    * </table>
    *
    * </center>
    *
    */
    template<typename Int, typename Dbl>
    class TestProjections {
      public:

        /**
         * `master` is the original lattice for which we want to calculate the
         * \f$M\f$ merit as defined above. `lattice` is the working lattice used for
         * intermediate calculations; all the projections from `master` are stored in
         * `lattice` before calling the test. `test` is the lattice test to be
         * applied on the different projections. `d` is the last element of array
         * `td`, which gives the maximal dimensions for the different projections.
         * `td[0]` and `td[1]` gives the minimal and the maximal dimensions for the
         * test applied on successive dimensions. For example, if `d` = 3 and `td` =
         * [2, 32, 16, 12], then the test will be applied on all successive
         * dimensions \f$2 \le t \le32\f$, on all possible 2-dimensional projections
         * for dimensions \f$\le16\f$, and on all possible 3-dimensional projections
         * for dimensions \f$\le12\f$. The results of the tests will be outputted on
         * `Writer`.
         */
        TestProjections (LatticeTester::IntLattice<Int, Int, Dbl, Dbl>
            *master, LatticeTester::IntLattice<Int, Int, Dbl, Dbl> 
            *lattice, LatticeTest<Int, Dbl> *test, int td[], int d);

        /**
         * Destructor.
         */
        ~TestProjections ();

        /**
         * Sends the output to `rw`. If this method is not called, output will
         * be sent to standard output.
         */
        void setOutput (LatticeTester::Writer<Int> * rw);

        /**
         * If `flag` is `true`, the tests will be applied on the *dual*
         * lattice; if it is `false`, the tests will be applied on the *primal*
         * lattice.
         */
        void setDualFlag (bool flag);

        /**
         * Sets the value of the <tt>m_invertF</tt> flag. If `invertF` is
         * `true`, the inverse of the length of the shortest vector will be
         * printed in the results. Otherwise, the length itself will be
         * printed.
         */
        void setInvertFlag (bool flag);

        /**
         * If flag is `true`, the value of the merit will be printed for all
         * projections, otherwise not. This flag should be set `false` in the
         * case of the `seek*` programs, and `true` otherwise.
         */
        void setPrintF (bool flag);

        /**
         * Builds the basis (and dual basis) of the projection `proj` for this
         * lattice. The result is placed in the work lattice. The basis is
         * triangularized to form a proper basis. For example, if \f$d=2\f$ and
         * `indices` = \f$[1, 4]\f$, then a 2-dimensional basis is built using
         * coordinates 1 and 4 of the master basis.
         */
        void build (const LatticeTester::Coordinates & proj);

        /**
         * Calculates the \f$M_{t_1, …, t_d}\f$ merit by running all the tests
         * for the different projections. If set `true`, `stationary` means
         * that the lattice is known to be *dimension stationary*. If set
         * `true`, the flag `last` means that only the projections which
         * include the last dimension of the lattice will be considered. This
         * is good for performing incremental searches for good lattices.
         * `minVal` is the minimal value of the normalized merit in each
         * dimension for a lattice to be considered. This method returns the
         * worst value of the merit over all projections.
         */
        double run (bool stationary, bool last, double minVal[]);

        /**
         * As method `run` above, but with the weights `weights`.
         */
        double run (bool stationary, bool last, double minVal[],
            const LatticeTester::Weights & weights);

        /**
         * Calculates the number of projections, given the parameters of this
         * object. The flag `stationary` must be set `true` if the lattice is
         * *dimension stationary*. If set `true`, the flag `last` means that
         * only the projections which include the last dimension of the lattice
         * will be considered.
         */
        int calcNumProjections (bool stationary, bool last);

        /**
         * Returns the number of projections considered for the last call to
         * `run`. If the test was interrupted because the lattice was rejected
         * as uninteresting, the number returned will be less than the total
         * number of projections.
         */
        int getNumProjections () { return m_numproj; }

        /**
         * Returns a pointer to the MeritProj pointed by m_meritproj.
         * */
        MeritProj * getMeritProj () { return m_meritproj; }

      protected:

        /**
         * Run only for projections given by the iterator `projit`.
         */
        double run (ProjIterator & projit, double minVal[],
            const LatticeTester::Weights & weights);

        /**
         * Lattice on which the \f$M_{t_1, …, t_d}\f$ merit will be calculated.
         */
        LatticeTester::IntLattice<Int, Int, Dbl, Dbl>* m_master;

        /**
         * Working lattice which is used to test the different projections.
         */
        LatticeTester::IntLattice<Int, Int, Dbl, Dbl>* m_lattice;

        /**
         * The lattice test used to calculate the merit.
         */
        LatticeTest<Int, Dbl>* m_test;

        /**
         * Weights.
         */
        double* m_weightsTemp;

        /**
         * If `true`, the test is applied on the dual lattice, otherwise on the
         * primal lattice.
         */
        bool m_dualF;

        /**
         * If `true`, the inverse length of the shortest vector is printed for
         * all projections, otherwise the length is printed.
         */
        bool m_invertF;

        /**
         * If `true`, the value of the merits is printed for all projections,
         * otherwise not.
         */
        bool m_printF;

        /**
         * If `true`, the square root of the values of the merit is printed
         * after the test; otherwise the values themselves are printed.
         */
        bool m_racF;

        /**
         * The maximal dimensions for each kind of projections.
         * <tt>m_td[0]</tt> is the minimal dimension for successive dimensions.
         * <tt>m_td[1]</tt> is the maximal dimension for successive dimensions
         * (1-dimensional projections), <tt>m_td[2]</tt> is the maximal
         * dimension for 2-dimensional projections, <tt>m_td[3]</tt> is the
         * maximal dimension for 3-dimensional projections, and so on. The last
         * element, <tt>m_td[m_d]</tt>, is the maximal dimension for
         * <tt>m_d</tt>-dimensional projections.
         */
        int *m_td;

        /**
         * The number of kinds of projections. Also the number of elements of
         * array <tt>m_td</tt> is <tt>m_d + 1</tt>.
         */
        int m_d;

        /**
         * The number of projections.
         */
        int m_numproj;

      private:

        /**
         * Output will be written on <tt>m_writer</tt>. By default, it is written on
         * standard output.
         */
        LatticeTester::Writer<Int> *m_writer;
        bool m_wrFlag;

        /**
         * Stores the figures of merit after `run` is executed.
         * */
        MeritProj* m_meritproj;
    };

  //===========================================================================

  template<typename Int, typename Dbl>
    TestProjections<Int, Dbl>::TestProjections (LatticeTester::IntLattice<Int, Int, 
        Dbl, Dbl> * master, LatticeTester::IntLattice<Int, Int, Dbl,
        Dbl> * lattice, LatticeTest<Int, Dbl> * test, int td[], int d)
    {
      if (d <= 0)
        LatticeTester::MyExit(1, "   TestProjections:   d <= 0");
      for (int i = 2; i <= d; i++) {
        if (td[i] >= 100)
          LatticeTester::MyExit(1, "   TestProjections:   td[i] >= 100");
      }
      m_d = d;
      m_td = new int[d + 1];
      memcpy (m_td, td, (d + 1) * sizeof (int));
      int maxDim = 0;
      for (int i = 0; i <= d; i++)
        maxDim = std::max (maxDim, m_td[i]);
      m_weightsTemp = new double[maxDim + 1];
      for (int i = 0; i <= maxDim; i++)
        m_weightsTemp[i] = 1.0;
      m_master = master;
      m_lattice = lattice;
      m_test = test;
      m_dualF = test->getDualFlag ();
      m_printF = true;
      m_invertF = test->getInvertFlag ();
      if ((LatticeTester::SPECTRAL == test->getCriterion() && LatticeTester::L2NORM == master->getNorm ())
          || (LatticeTester::BEYER == test->getCriterion()))
        m_racF = true;
      else
        m_racF = false;
      m_writer = new LatticeTester::WriterRes<Int> (&std::cout);
      m_wrFlag = true;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    TestProjections<Int, Dbl>::~TestProjections ()
    {
      delete[] m_td;
      delete[] m_weightsTemp;
      if (m_wrFlag)
        delete m_writer;
      delete m_meritproj;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void TestProjections<Int, Dbl>::setOutput (LatticeTester::Writer<Int> * rw)
    {
      delete m_writer;
      m_wrFlag = false;
      m_writer = rw;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void TestProjections<Int, Dbl>::setDualFlag (bool dual)
    {
      m_dualF = dual;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void TestProjections<Int, Dbl>::setPrintF (bool flag)
    {
      m_printF = flag;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void TestProjections<Int, Dbl>::build (const LatticeTester::Coordinates & proj)
    {
      m_master->buildProjection (m_lattice, proj);
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    int TestProjections<Int, Dbl>::calcNumProjections (bool stationary, bool forceLast)
    {
      // Les projections sur dimensions successives
      m_numproj = m_td[1] - m_td[0] + 1;

      // Les projections sur des dimensions successives qui ne contiennent pas 0
      //comme indice
      if (!stationary) {
        for (int order = 2; order <= m_d; order++) {
          int maxCoord = std::min (m_td[order], m_master->getDim ());
          if (maxCoord <= order)
            continue;
          for (ProjIteratorSuccCoords projit(2, maxCoord, order, order, stationary, forceLast); projit; ++projit)
            m_numproj++;
        }
      }

      // test for non-successive coordinates
      // loop over projection orders
      for (int order = 2; order <= m_d; order++) {
        int maxCoord = std::min (m_td[order], m_master->getDim ());
        if (maxCoord <= order)
          continue;
        for (ProjIteratorNonSuccCoords projit(1, maxCoord, order, order, stationary, forceLast); projit; ++projit)
          m_numproj++;
      }

      //! if (false == stationary) {
      //!    // Les autres projections sur dimensions successives
      //!    for (int i = 2; i <= m_d; i++) {
      //!       int dim = std::min (m_td[i], m_master->getMaxDim ());
      //!       initSuccDims (i, dim, forceLast);
      //!       do {
      //!          m_numproj++;
      //!       } while (nextSuccDims (stationary, forceLast, i, dim));
      //!    }
      //! }

      //! // Les projections sur dimensions non successives
      //! for (int i = 2; i <= m_d; i++) {
      //!    int dim = std::min (m_td[i], m_master->getMaxDim ());
      //!    if (initNonSuccDims (i, dim, forceLast)) {
      //!       do {
      //!          m_numproj++;
      //!       } while (nextNonSuccDims (stationary, forceLast, i, dim));
      //!    }
      //! }
      return m_numproj;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    double TestProjections<Int, Dbl>::run (bool stationary, bool forceLast, double minVal[])
    {
      return run(stationary, forceLast, minVal, unitWeights);
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    double TestProjections<Int, Dbl>::run (bool stationary, bool forceLast, double minVal[],
        const LatticeTester::Weights& weights)
    {
      m_meritproj = new MeritProj(calcNumProjections(stationary, forceLast));
      // we assume that m_td[1] >= m_td[i] for all i
      int maxDim = m_td[1];
      double merit = 1.0e100;

      //! std::cout << "minVal: " << minVal[0];
      //! for (int j = 1; j < maxDim; j++)
      //!    std::cout << "," << minVal[j];
      //! std::cout << std::endl;

      // Le test pour les projections sur dimensions successives
      m_lattice->buildBasis (m_td[0] - 1);
      int minDim = m_td[0];
      m_test->setDualFlag (m_dualF);
      m_test->setMaxAllDimFlag (true);

      // set the temporary weights for successive dimensions
      for (ProjIteratorSuccCoords projit(0, maxDim-1, 0, maxDim-1, true, false);
          projit; ++projit) {
        if ((int)projit->size() >= minDim)
          m_weightsTemp[projit->size()] = weights.getWeight(*projit);
        else
          m_weightsTemp[projit->size()] = 1;
      }

      m_test->test (minDim, maxDim, minVal, m_weightsTemp);
      // ATTENTION: si le test s'est terminé prématurément parce que le réseau
      // est mauvais, les valeurs de mérites ci-après sont n'importe quoi.
      merit = m_test->getMerit().getST (minDim, maxDim);
      m_numproj = maxDim - minDim + 1;
      for (int i = 0; i < m_numproj; i++) {
        (*m_meritproj)[i] = m_test->getMerit()[minDim+i];
        m_meritproj->getMerit(i) = m_test->getMerit().getMerit(minDim+i);
        std::string a("{");
        for (int j = 0; j < minDim; j++) {
          if (j != 0) a += ",";
          a += std::to_string(j);
        }
        for (int j = 0; j < i; j++) a += "," + std::to_string(j+minDim);
        a += "}";
        m_meritproj->setCoord(i, a);
      } 

      if (m_printF) {
        //std::cout << "------------------------------------------" << std::endl;
        //std::cout << " a = " << m_lattice->toStringCoef () << std::endl;
        m_writer->writeString (" ");
        m_writer->getStream() << weights;
        m_writer->newLine ();
        m_writer->writeString (
            m_test->getMerit ().toString(minDim, maxDim, m_racF, m_invertF));
        printMerit (m_racF, merit, m_writer);
      }
      if (merit < minVal[maxDim] || m_d <= 1)
        return merit;

      /*
         if (!stationary) {
      // test for successive coordinates
      // loop over projection orders
      for (int order = 2; order <= m_d; order++) {

      int maxCoord = min (m_td[order], m_master->getDim ());
      // if maxCoord < order, there are no projections to consider
      // if maxCoord == order, the merit has already been computed above
      if (maxCoord <= order)
      continue;

      ProjIteratorSuccCoords projit(1, maxCoord-1, order, order, stationary,
      forceLast);
      merit = std::min(merit, run(projit, minVal, weights));

      if (merit < minVal[order])
      return merit;

      if (m_printF)
      printMerit (m_racF, merit, m_writer);
      }
      }
      */

      // Testing projections of size order
      for (int order = 2; order <= m_d; order++) {

        int maxCoord = std::min (m_td[order], m_master->getDim ());
        // if maxCoord <= order, there are no projections with non-successive
        // indices to consider
        if (maxCoord <= order)
          continue;

        // Iterates over sets of `order` non successive coordinates between `0` 
        // and `maxCoord - 1`
        ProjIteratorNonSuccCoords projit(0, maxCoord-1, order, order, stationary,
            forceLast);
        merit = std::min(merit, run(projit, minVal, weights));

        // If the lattice is not dimension stationary, we look for successive
        // coordinates projections too
        if (!stationary) {
          ProjIteratorSuccCoords projit(1, maxCoord-1, order, order, stationary,
              forceLast);
          merit = std::min(merit, run(projit, minVal, weights));
        }

        if (merit < minVal[order])
          return merit;

        if (m_printF)
          printMerit (m_racF, merit, m_writer);
      }

      m_test->getMerit ().setWorstMerit (merit);
      return merit;
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    double TestProjections<Int, Dbl>::run (ProjIterator& projit, double minVal[], 
        const LatticeTester::Weights& weights)
    {
      std::ostringstream os;
      double minMerit = 1.0e100;
      while (projit) {
        int dim = (int)projit->size();
        int maxCoord = std::min (m_td[dim], m_master->getDim ());
        if (maxCoord <= dim) // if equal, already computed
          continue;

        if (m_printF)
          m_writer->writeString (formatIndices(*projit));

        m_numproj++;

        m_weightsTemp[dim] = weights.getWeight(*projit);

        m_master->buildProjection (m_lattice, *projit);

        m_test->test (dim, dim, minVal, m_weightsTemp);

        double len = m_test->getMerit ().getMerit (dim);

        double curMerit = m_test->getMerit ().getNormVal (dim);
        if (m_printF)
          printMerLen (m_invertF, m_racF, len, curMerit, m_writer);
        if (curMerit < minVal[dim])
          return curMerit;

        os.str(std::string());
        os << *projit;
        (*m_meritproj)[m_numproj-1] = m_test->getMerit()[dim];
        m_meritproj->getMerit(m_numproj-1) = m_test->getMerit().getMerit(dim);
        m_meritproj->getCoord(m_numproj-1) = os.str();

        minMerit = std::min (minMerit, curMerit);

        // incrementing stuff
        ++projit;
      }
      return minMerit;
    }

}
#endif

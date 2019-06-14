#ifndef LATMRG_FIGUREOFMERIT_H
#define LATMRG_FIGUREOFMERIT_H

#include <list>

#include "latmrg/Projections.h"

namespace LatMRG {

  /**
   * This class is intended to store everything you might need after performing
   * computations on a lattice. It is intended as a return type for functions
   * actually performing the computations.
   * */
  template<typename Lat>
    class FigureOfMerit {
      typedef typename Lat::Float Dbl;

      private:
      /**
       * The lattice that has been tested.
       * */
      Lat m_lattice;

      /**
       * `true` if this contains a merit value for all projections, `false`
       * otherwise.
       * */
      bool m_finished;

      /**
       * The projections used in this object.
       * */
      Projections m_projSet;

      /**
       * The number of projections this has been tested over.
       * */
      long m_numProj;

      /**
       * The merit for each test in the same sequence as the projections in
       * `m_projSet`.
       * */
      std::vector<Dbl> m_merits;

      /**
       * The shortest vectors for each projection in the same sequence as in
       * `m_projSet`.
       * */
      std::vector<typename Lat::IntVec> m_vectors;

      /**
       * The merit of the lattice.
       * */
      Dbl m_merit;

      /**
       * Indice of the best or worst value in `m_merits`.
       * */
      long m_best_worst;

      public:

      FigureOfMerit() {
        finished = false;
      }

      /**
       * Initializes the object for `lattice` and `proj`.
       * This won't contain any merit value.
       * */
      FigureOfMerit(const Lat& lattice, const Projections& proj):
        m_lattice(lattice), m_projSet(proj) {
          m_finished = false;
        }

      /**
       * Computes the minimum, if `strat == "min"`, or the maximum, if
       * `strat == "max"`, of `m_merits` and stores it in `m_merit`.
       * */
      void computeMerit(const std::string strat) {
        if (!m_finished) {
          m_merit = -1;
          m_best_worst = 0;
          return;
        }
        if (strat == "min") {
          m_merit = m_merits[0];
          m_best_worst = 0;
          for (unsigned int i = 1; i < m_merits.size(); i++) {
            if (m_merit > m_merits[i]) {
              m_merit = m_merits[i];
              m_best_worst = i;
            }
          }
        } else if (strat == "max") {
          m_merit = m_merits[0];
          m_best_worst = 0;
          for (unsigned int i = 1; i < m_merits.size(); i++) {
            if (m_merit < m_merits[i]) {
              m_merit = m_merits[i];
              m_best_worst = i;
            }
          }
        }
      }

      /**
       * Sets `m_merits` to `merits` and `m_vectors` to `vectors`.
       * */
      void addMerit(std::vector<Dbl>& merits, std::vector<typename Lat::IntVec>& vectors) {
        m_merits = merits;
        m_vectors = vectors;
      }

      /**
       * Sets finished to `true`.
       * */
      void setFinished() { m_finished = true;}

      /**
       * Returns the shortest vector in the lattice for the worst projection.
       * */
      typename Lat::IntVec worstVect() {return m_vectors[m_best_worst];}

      /**
       * */
      LatticeTester::Coordinates worstProj(){
        m_projSet.resetDim();
        for (int i = 0; i < m_best_worst; i++) ++m_projSet;
        return m_projSet.getProj();
      }


      /**
       * Returns `m_finished`.
       * */
      bool finished() { return m_finished;}

      /**
       * Returns `m_merit`.
       * */
      Dbl getMerit(){
        return m_merit;
      }

      /// Returns the string associated with this test.
      std::string getLattice() {
        return m_lattice.toString();}

      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator==(const FigureOfMerit<Lat>& right) {
        return m_merit == right.m_merit;
      }
      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator< (const FigureOfMerit<Lat>& right) {
        return m_merit < right.m_merit;
      }
      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator> (const FigureOfMerit<Lat>& right) {
        return right.m_merit < m_merit;
      }
      /**
       * Compares `m_merit`. Uses `operator>`.
       * */
      inline bool operator<=(const FigureOfMerit<Lat>& right) {
        return !operator>(right);
      }
      /**
       * Compares `m_merit`. Uses `operator<`.
       * */
      inline bool operator>=(const FigureOfMerit<Lat>& right) {
        return !operator<(right);
      }
      /**
       * Compares `m_merit`. Uses `operator==`.
       * */
      inline bool operator!=(const FigureOfMerit<Lat>& right) {
        return !operator==(right);
      }

    };

  /**
   * This class stores a certain number of `FigureOfMerit` objects.
   * New elements are only added if the list is not full or if they have a better
   * merit than what is on the list.
   * */
  template<typename Lat>
    class MeritList{
      typedef typename Lat::Float Dbl;
      public:
      /**
       * Creates a `MeritList` that will hold at most `maxLength` `FigureOfMerit`.
       * */
      MeritList(unsigned long maxLength, bool best){
        m_max = maxLength;
        m_merit = Dbl(0);
        m_best = best;
      }

      /**
       * Adds `test` in `m_tests` if there is space left or if `test` as a
       * higher merit.
       * */
      void add(FigureOfMerit<Lat> test) {
        if (!test.finished()) return;
        if ((m_tests.size() >= m_max) && !(m_best ^ (test > m_tests.back()))) {
          m_tests.pop_back();
        } else if (m_tests.size() >= m_max) return;
        posInsert(test);
      }

      unsigned long getSize(){
        return m_tests.size();
      }

      std::list<FigureOfMerit<Lat>>& getList() {return m_tests;}

      Dbl getMerit(){
        return m_merit;
      }

      private:

      /**
       * The tests stored in this object.
       * */
      std::list<FigureOfMerit<Lat>> m_tests;

      /**
       * The maximal number of tests this object will store.
       * */
      unsigned long m_max;

      /**
       * The minimum merit stored in this object.
       * */
      Dbl m_merit;

      /**
       * Indicates wether or not this object keeps best or worst lattices.
       * */
      bool m_best;

      /**
       * Finds the right position in the vector and inserts the new test.
       * */
      void posInsert(FigureOfMerit<Lat>& test){
        if (m_tests.size() == 0) {
          m_tests.push_back(test);
          return;
        }
        for (auto rit = m_tests.crbegin(); rit != m_tests.crend(); rit++) {
          if (m_best) {
            if ((test < *rit)) {
              m_tests.insert(rit.base(), test);
              if (rit == m_tests.crbegin() && m_tests.size() == m_max) m_merit = test.getMerit();
              return;
            }
          } else {
            if ((test > *rit)) {
              m_tests.insert(rit.base(), test);
              if (rit == m_tests.crbegin() && m_tests.size() == m_max) m_merit = test.getMerit();
              return;
            }
          }
        }
        m_tests.insert(m_tests.cbegin(), test);
      }
    };

} // end namespace LatMRG
#endif

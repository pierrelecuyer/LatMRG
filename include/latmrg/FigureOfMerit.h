#ifndef LATMRG_FIGUREOFMERIT_H
#define LATMRG_FIGUREOFMERIT_H

namespace LatMRG {

  /**
   * This class is intended to store everything you might need after performing
   * computations on a lattice. It is intended as a return type for functions
   * actually performing the computations.
   *
   * \todo add stuff here.
   * This should contain:
   * - It should be a template class, the lattice type behing what is the template
   *   parameter
   * - The projections set
   * - The merit on each projection
   * - The shortest vector in each projection
   * - A recipe for merit computation
   * - a method to compute merits
   * - A set of weights to perform a figure of merit computation
   * - ? a string describing the test
   * */
  class FigureOfMerit {
    public:
      Test(const std::string& lattice, const DblVec& results) {
        m_lattice = lattice;
        m_results = DblVec(results);
        computeMerit();
      }

      /**
       * Computes the minimum of `m_results` and stores it in `m_merit`.
       * */
      void computeMerit() {
        m_merit = m_results[0];
        for (int i = 1; i < m_results.length(); i++) {
          if (m_merit > m_results[i]) m_merit = m_results[i];
        }
      }

      /**
       * Returns `m_merit`.
       * */
      Dbl getMerit(){
        return m_merit;
      }

      /// Returns the string associated with this test.
      std::string getLattice() {return m_lattice;}

      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator==(const Test& right) {
        return m_merit == right.m_merit;
      }
      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator< (const Test& right) {
        return m_merit < right.m_merit;
      }
      /**
       * Compares `m_merit`. "Implements" the comparison.
       * */
      inline bool operator> (const Test& right) {
        return right.m_merit < m_merit;
      }
      /**
       * Compares `m_merit`. Uses `operator>`.
       * */
      inline bool operator<=(const Test& right) {
        return !operator>(right);
      }
      /**
       * Compares `m_merit`. Uses `operator<`.
       * */
      inline bool operator>=(const Test& right) {
        return !operator<(right);
      }
      /**
       * Compares `m_merit`. Uses `operator==`.
       * */
      inline bool operator!=(const Test& right) {
        return !operator==(right);
      }


    private:
      /**
       * The lattice that has been tested.
       * */
      std::string m_lattice;

      /**
       * The normalized results of the test.
       * */
      DblVec m_results;

      /**
       * The merit of the lattice.
       * */
      Dbl m_merit;
  };

  /**
   * This class stores a certain number of `FigureOfMerit` objects.
   * New elements are only added if the list is not full or if they have a better
   * merit than what is on the list.
   * */
  class MeritList{
    public:
      /**
       * Creates a `TestList` that will hold at most `maxLength` `Test`.
       * */
      TestList(unsigned long maxLength, bool best){
        m_max = maxLength;
        m_merit = Dbl(0);
        m_best = best;
      }

      /**
       * Adds `test` in `m_tests` if there is space left or if `test` as a
       * higher merit.
       * */
      void add(Test& test) {
        if ((m_tests.size() >= m_max) && !(m_best ^ (test > m_tests.back()))) {
          m_tests.pop_back();
        } else if (m_tests.size() >= m_max) return;
        posInsert(test);
      }

      unsigned long getSize(){
        return m_tests.size();
      }

      std::list<Test> getList() {return m_tests;}

      Dbl getMerit(){
        return m_merit;
      }

    private:

      /**
       * The tests stored in this object.
       * */
      std::list<Test> m_tests;

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
      void posInsert(Test& test){
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

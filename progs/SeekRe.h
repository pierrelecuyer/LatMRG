#ifndef SEEKRE_H
#define SEEKRE_H

/**
 * This file contains functions, classes and code used in SeekRe. This file is
 * used to remove some of the clutter in SeekRe.cc. This is because some of the
 * functions needed already exist but need to be rewritten.
 * */
#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/ParamReaderExt.h"

typedef NTL::ZZ Int;
typedef double Dbl;
typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;
typedef NTL::vector<Dbl> DblVec;

using LatticeTester::IntLattice;

namespace LatMRG {

  /**
   * This class stores a lattice with the results to a test. This class is used
   * to agregate the information of a test with its lattice. This class computes
   * the merit of the lattice as the minimum of the results table given to it
   * and implements comparison operators to compare the merit of multiple
   * lattices.
   * */
  class Test{
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
   * This class stores a certain number of `Test` objects. When adding a `Test`
   * to this class, it first checks
   * */
  class TestList{
    public:
      /**
       * Creates a `TestList` that will hold at most `maxLength` `Test`.
       * */
      TestList(unsigned long maxLength){
        m_max = maxLength;
      }

      /**
       * Adds `test` in `m_tests` if there is space left or if `test` as a
       * higher merit.
       * */
      void add(Test& test) {
        if ((m_tests.size() >= m_max) && (test > m_tests.back())) {
          m_tests.pop_back();
        } else if (m_tests.size() >= m_max) return;
        posInsert(test);
      }

      unsigned long getSize(){
        return m_tests.size();
      }

      std::list<Test> getList() {return m_tests;}

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
       * Finds the right position in the vector and inserts the new test.
       * */
      void posInsert(Test& test){
        if (m_tests.size() == 0) {
          m_tests.push_back(test);
          return;
        }
        for (auto rit = m_tests.crbegin(); rit != m_tests.crend(); rit++) {
          if (test < *rit) {
            m_tests.insert(rit.base(), test);
            return;
          }
        }
        m_tests.insert(m_tests.cbegin(), test);
      }
  };

  /*
   * A more flexible projection class. This class contains all the information
   * of the projections to test, this includes the dimensions for each type of
   * projection. 
   *
   * This only builds projections with the coordinate 0 because the projections
   * of a MRG lattice are all equivalent to a projection with 0 as a coordinate
   * */
  class Projections {
    public:
      Projections(int dimProj, int minDim, std::vector<std::size_t>& projDim){
        std::cout << "dimProj: " << dimProj;
        std::cout << "\nminDim: " << minDim << "\n";
        for (auto it = projDim.begin(); it != projDim.end(); it++)
          std::cout << *it << " ";
        std::cout << "\n";
        m_numDim = dimProj;
        m_projDim.resize(dimProj);
        m_minDim = minDim;
        for (int i = 0; i < dimProj; i++) {
          m_projDim[i] = projDim[i];
        }
        m_currentDim = 1;
      }

      LatticeTester::Coordinates next() {
        if (m_curProj.empty()) {
          for (std::size_t i = 0; i<(unsigned)m_minDim; i++) m_curProj.push_back(i);
        }
        return LatticeTester::Coordinates(m_curProj);
      }

      LatticeTester::Coordinates getProj() {
        return LatticeTester::Coordinates(m_curProj);
      }

      /**
       * Returns true if the current projection is the last one.
       * */
      bool end() {
        if (m_curProj.empty()) return false;
        return true;
      }

      void resetDim(int dim = 0) {
        m_currentDim = dim;
        m_curProj.clear();
      }

      int getDim() { return m_currentDim;}

    private:
      int m_numDim;

      std::vector<std::size_t> m_curProj;

      int m_currentDim;

      int m_minDim;

      std::vector<std::size_t> m_projDim;
  };

}

#endif

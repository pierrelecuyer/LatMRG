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
   * A projection class that does exactly what we need when generating
   * projections when testing a generator.
   *
   * This is initialized containing nothing and can iteratively return
   * projections by calling `next()`. This first goes a through a list of
   * sequential projections and then gives a list of projections with non
   * sequential indices.
   * */
  class Projections {
    private:
      /*
       * The number of dimension this has projections on. This can be misleading
       * because what we consider to be the first dimension is in fact
       * projections on dimensions with sequential indices.
       * */
      int m_numDim;

      /*
       * The current projection stored in this object.
       * */
      std::vector<std::size_t> m_curProj;

      /*
       * The dimension of the current projection. If this is 1, this does not
       * represent the actual dimension of the projection, but only indicates
       * that the indices are sequential. This is not what is returned by getDim().
       * */
      int m_currentDim;

      /*
       * Every projection will include at least one coordinate equal or greater
       * than `m_minDim-1`. This helps when studying full period generators.
       * */
      int m_minDim;

      /*
       * `m_projDim[i]` is the maximal indice used by projections when
       * `m_currentDim == i+1`.
       * */
      std::vector<std::size_t> m_projDim;

    public:

      /**
       * Creates a Projections set with sequential projections and non sequential
       * projections in dimensions 2 through `dimProj`. Projections in dimension
       * `i` use indices up to `projDim[i-1]` and sequential projections use
       * indices up to `projDim[0]`
       * */
      Projections(int dimProj, int minDim, std::vector<std::size_t>& projDim){
        m_numDim = dimProj;
        m_projDim.resize(dimProj);
        m_minDim = minDim;
        for (int i = 0; i < dimProj; i++) {
          m_projDim[i] = projDim[i];
        }
        resetDim();
      }

      /**
       * If `dim==0` this will return true if this calling next will not get a
       * new projection.
       * If `dim==1` this will return true if the current projection dimension
       * is over.
       * */
      bool end(int dim = 0) {
        if (m_currentDim == 0 || m_curProj.size() == 0) return false;
        //std::cout << LatticeTester::Coordinates(m_curProj) << "\n";
        bool cond = m_curProj.back() == m_projDim[m_currentDim-1];
        //std::cout << "cond: " << cond << " dim: " << dim << "\n";
        // Exhausted all options
        if (dim == 0 && m_currentDim == m_numDim) {
          if (m_numDim == 1 || m_numDim == 2) return cond;
          // Second coordinate is in last possible spot (see how we increment in
          // `next()`)
          //std::cout << "dim 0: " << (m_curProj[1] == (m_projDim[m_numDim-1]-m_numDim+2)) << "\n";
          return m_curProj[1] == (m_projDim[m_numDim-1]-m_numDim+2);
        }
        if (dim == 1) {
          if (m_currentDim == 1 || m_currentDim == 2) return cond;
          //std::cout << "dim 1: " << (m_curProj[1] == (m_projDim[m_currentDim-1]-m_currentDim+2)) << "\n";
          return m_curProj[1] == (m_projDim[m_currentDim-1]-m_currentDim+2);
        }
        return false;
      }
      
      /**
       * This returns the next projection for this set of projections.
       *
       * If this object is not initialized, this fills it with the first
       * projection with indices from `0` to `m_minDim-1`. It then adds
       * coordinates up to `m_projDim[0]`. It then returns coordinates 
       * */
      LatticeTester::Coordinates next() {
        if (end()) return LatticeTester::Coordinates(m_curProj);
        // Filling the vector for the first time, assuming it is empty.
        if (m_currentDim == 0) {
          for (std::size_t i = 0; i<(unsigned)m_minDim; i++) m_curProj.push_back(i);
          m_currentDim++;
          return LatticeTester::Coordinates(m_curProj);
        }
        // Checking dimension change
        if (end(1)) {
          m_currentDim++;
          m_curProj.resize(0);
        }
        // Dimension 1: add a coordinate to the vector
        if (m_currentDim == 1) {
          m_curProj.push_back(m_curProj.size());
        } else if (m_curProj.size() == 0) {
          // Generating the first vector for a non-sequential projection in dim
          // m_currentDim
          for (std::size_t i = 0; i<(unsigned)(m_currentDim-1); i++)
            m_curProj.push_back(i);
          // Projections without a coordinate greater or equal to m_minDim-1 are
          // irrelevant as are fully sequential projections
          m_curProj.push_back((m_currentDim>m_minDim-1)?
              (unsigned)m_currentDim:(unsigned)(m_minDim-1));
        } else if (m_currentDim == 2) {
          // We know we can safely increment the last coordinate
          m_curProj[1] = m_curProj[1] + 1;
        } else  {
          // Incrementing a vector that can be incremented
          int i = m_currentDim-1;
          while (m_curProj[i] == m_projDim[m_currentDim-1]-m_currentDim+1+i) i--;
          m_curProj[i]++;
          for(i+=1; i<m_currentDim; i++) {
            m_curProj[i] = m_curProj[i-1]+1;
          }
          if (m_curProj.back() < (unsigned)(m_minDim-1)) m_curProj.back() = m_minDim-1;
        }
        return LatticeTester::Coordinates(m_curProj);
      }

      LatticeTester::Coordinates getProj() {
        return LatticeTester::Coordinates(m_curProj);
      }

      void resetDim(int dim = 0) {
        if (dim > m_numDim) return;
        m_currentDim = dim;
        m_curProj.resize(0);
      }

      int getDim() { return m_curProj.size();}

  };

}

#endif

#ifndef LATMRG_PROJECTIONS_H
#define LATMRG_PROJECTIONS_H

#include <vector>

#include "latticetester/Coordinates.h"

namespace LatMRG {

  /**
   * This class permits one to generate the set of projections that we need
   * to examine when computing a figure of merit.
   * After initializing the `Projections` object, the desired successive projections are generated
   * and returned iteratively by calling `next()`.
   *
   * FROM PIERRE: Unclear why we need this class.  It seems that `CoordinateSets` is
   * doing the same thing already.  What is the difference?
   * Also, the design is not very natural, not clearly explained, and can lead to errors.
   *
   *
   * In the current implementation, this first goes a through a list of
   * projections over sequential coordinates starting at 0, and then
   * through a list of projections over non-successive coordinates.
   *
   * FROM PIERRE: This order is too much of a special case, too rigid and probably not optimal!
   * In practice, it may be more efficient to look at the pairs, triples, etc., in this order.
   * One should have more flexibility for the choice of ordering!
   * Also, the current design of this class is too complicated, not very elegant!
   *
   * When iterating through this object, it returns `LatticeTester::Coordinates`
   * objects. These contain coordinates starting from `0` to be used when
   * accessing elements in vectors. Take care, some of the functions in this
   * object take integers starting at 1 that are transformed by the object.
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
       * This is also the smallest dimension that will be tested for projections
       * on sequential coordinates.
       * */
      int m_minDim;

      /*
       * `m_projDim[i]` is the maximal indice used by projections when
       * `m_currentDim == i+1`.
       * */
      std::vector<int> m_projDim;

    public:

      /**
       * Creates a Projections set. It is necessary that `dimProj == length(projDim)`.
       *
       * This set will contain projections with sequential indices in dimensions
       * `minDim` to `projDim[0]+1`. That means the projections will include
       * sets `{0, ..., minDim-1}` through `{0, ..., projDim[0]}`.
       *
       * When the length of `projDim` is greater than 1, this set will contain
       * non sequential projections in dimensions 2 through `dimProj`.
       * Projections in dimension `i` use `i` indices in the range
       * `[0, projDim[i-1]]` with at least one of them at least `minDim-1`.
       *
       * FROM PIERRE: It seems that `projDim` contains the values of t_1, t_2, ... t_d.
       * */
      Projections(int dimProj, int minDim, std::vector<int>& projDim);

      /**
       * If `dim==0` this will return true if this calling next will not get a
       * new projection.
       * If `dim==1` this will return true if the current projection dimension
       * is over.
       * */
      bool end(int dim = 0);

      /**
       * This returns the next projection for this set of projections.
       *
       * If this object is not initialized, this fills it with the first
       * projection with indices from `0` to `m_minDim-1`. It then adds
       * coordinates up to `m_projDim[0]`. It then returns coordinates
       * */
      LatticeTester::Coordinates next();

      /**
       * Returns the current projection the object is at.
       * */
      LatticeTester::Coordinates getProj();

      void resetDim(int dim = 0);

      /**
       * Returns the dimension of the current projections. This does not return
       * the value of `m_currentDim`
       * */
      int getDim();

      int minDim() {return m_minDim;}
      int maxDim() {return m_projDim[0]+1;}
      int numProj() {return m_numDim;}
      std::vector<int>& projDim() {return m_projDim;}

      /**
       * Returns a string describing this object of the form
       * ```
       * Sequential: minDim to m_projDim[0]+1
       * Projections Dimension                  i
       *             Coordinates m_projDim[i-1]+1
       * ```
       * */
      std::string toString();

      /**
       * Returns the result of next(). This is for convenience.
       * */
      LatticeTester::Coordinates operator++();

  }; // End class Projections


//============================================================================
// Implementation
  
  Projections::Projections(int dimProj, int minDim, std::vector<int>& projDim){
      m_numDim = dimProj;
      m_projDim.resize(dimProj);
      m_minDim = minDim;
      for (int i = 0; i < dimProj; i++) {
        m_projDim[i] = projDim[i];
      }
      resetDim();
    }

    //============================================================================

    bool Projections::end(int dim) {
      if (m_currentDim == 0 || m_curProj.size() == 0) {
        if (dim == 1 && !(m_projDim[m_currentDim-1] >= m_currentDim)) return true;
        return false;
      }
      //std::cout << LatticeTester::Coordinates(m_curProj) << "\n";
      bool cond = (int)m_curProj.back() == m_projDim[m_currentDim-1];
      //std::cout << "cond: " << cond << " dim: " << dim << "\n";
      // Exhausted all options
      if (dim == 0 && m_currentDim == m_numDim) {
        if (m_numDim == 1 || m_numDim == 2) return cond;
        // Second coordinate is in last possible spot (see how we increment in
        // `next()`)
        //std::cout << "dim 0: " << (m_curProj[1] == (m_projDim[m_numDim-1]-m_numDim+2)) << "\n";
        return (int)m_curProj[1] == (m_projDim[m_numDim-1]-m_numDim+2);
      }
      if (dim == 1) {
        if (!(m_projDim[m_currentDim-1] >= m_currentDim)) return true;
        if (m_currentDim == 1 || m_currentDim == 2) return cond;
        //std::cout << "dim 1: " << (m_curProj[1] == (m_projDim[m_currentDim-1]-m_currentDim+2)) << "\n";
        return (int)m_curProj[1] == (m_projDim[m_currentDim-1]-m_currentDim+2);
      }
      return false;
    }

    //============================================================================

    LatticeTester::Coordinates Projections::next() {
      if (end()) return LatticeTester::Coordinates(m_curProj);
      // Filling the vector for the first time, assuming it is empty.
      if (m_currentDim == 0) {
        for (std::size_t i = 0; i<(unsigned)m_minDim; i++) m_curProj.push_back(i);
        m_currentDim++;
        return LatticeTester::Coordinates(m_curProj);
      }
      // Checking dimension change
      while (end(1)) {
        m_currentDim++;
        m_curProj.resize(0);
      }
      // Dimension 1: add a coordinate to the vector
      if (m_currentDim == 1) {
        while(m_curProj.size() < (unsigned)(m_minDim-1)) m_curProj.push_back(m_curProj.size());
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
        while ((int)m_curProj[i] == m_projDim[m_currentDim-1]-m_currentDim+1+i) i--;
        m_curProj[i]++;
        for(i+=1; i<m_currentDim; i++) {
          m_curProj[i] = m_curProj[i-1]+1;
        }
        if (m_curProj.back() < (unsigned)(m_minDim-1)) m_curProj.back() = m_minDim-1;
      }
      return LatticeTester::Coordinates(m_curProj);
    }

    //============================================================================

    LatticeTester::Coordinates Projections::getProj() {
      return LatticeTester::Coordinates(m_curProj);
    }

    //============================================================================

    void Projections::resetDim(int dim) {
      if (dim > m_numDim) return;
      m_currentDim = dim;
      m_curProj.resize(0);
    }

    //============================================================================

    int Projections::getDim() { return m_curProj.size();}

    //============================================================================

    std::string Projections::toString() {
      std::string out("Sequential: ");
      out += std::to_string(m_minDim) + " to " + std::to_string(m_projDim[0]+1) + "\n";
      if (m_numDim > 1) {
        out += "Projections: ";
        out += "Dimension   ";
        for (int i = 2; i <= m_numDim; i++) {
          if (i >= 10) out += "   ";
          else out += "    ";
          out += std::to_string(i);
        }
        out += "\n";
        out += "             Coordinates ";
        for (int i = 1; i < m_numDim; i++) {
          if (m_projDim[i] >= 9){
            if (m_projDim[i] >= 99) out += "  ";
            else out += "   ";
          }
          else out += "    ";
          out += std::to_string(m_projDim[i]+1);
        }
        out += "\n";
      }
      return out;
    }

    //============================================================================

    LatticeTester::Coordinates Projections::operator++() {
      return next();
    }

  
} // End namespace LatMRG
#endif

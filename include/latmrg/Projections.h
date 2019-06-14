#ifndef LATMRG_PROJECTIONS_H
#define LATMRG_PROJECTIONS_H

#include <vector>

#include "latticetester/Coordinates.h"

namespace LatMRG {
  /*
   * A projection class that does exactly what we need when generating
   * projections when testing a generator.
   *
   * This is initialized containing nothing and can iteratively return
   * projections by calling `next()`. This first goes a through a list of
   * sequential projections and then gives a list of projections with non
   * sequential indices.
   *
   * When iterating through this object, it returns `LatticeTester::Coordinates`
   * objects. These contain coordinates starting from `0` to be used when
   * accessing elements in vectors. Take care, some of the functions in this
   * objects take integers starting at 1 that are transformed by the object.
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
      std::vector<std::size_t> m_projDim;

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
       * */
      Projections(int dimProj, int minDim, std::vector<std::size_t>& projDim);

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
      std::vector<std::size_t>& projDim() {return m_projDim;}

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

} // End namespace LatMRG
#endif

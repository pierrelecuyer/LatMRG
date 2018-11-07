#ifndef MERITPROJ_H
#define MERITPROJ_H
#include "latmrg/Merit.h"

#include <string>
#include <iostream>

namespace LatMRG {

  /**
   * A class storing the merit information when projections in multiple
   * dimensions are considered. This behaves slightly differently then Merit in
   * the sense that instead of storing the merit of a LatticeTest after a call
   * of `test`, it stores the merit for all the projections considered when
   * calling TestProjections.run(). Currently, this class does not segment the 
   * figures of merit depending on dimension because it is solely used for 
   * printing.
   * */
  class MeritProj: public Merit {
    public:
      /**
       * Empty constructor. Does nothing.
       * */
      MeritProj(){};

      /**
       * Initializes an object of this class to store the figures of merit for
       * `numproj` projections.
       * */
      MeritProj(int numproj);

      /**
       * Sets the name of coordinates of the the i-th projection to `name`.
       * */
      void setCoord(int i, std::string name){m_coord[i] = name;}

      /**
       * Sets the name of coordinates to an empty string.
       * */
      void setCoord();

      /**
       * Gets the name of coordinates of the the i-th projection.
       * */
      std::string & getCoord(int i){return m_coord[i];}

      /**
       * Gets the number of projections this object can contain.
       * */
      int getNumProj(){return m_numproj;}

      /**
       * Returns a string representing this object. This can be used to print
       * the figures of merit for all projection that have been stored in this
       * object. This has the same options `rac` and `invert` as Merit.toString.
       * */
      std::string toString (bool rac, bool invert);

    private:

      /**
       * The number of projections that this object contains. This should at 
       * all times be the same as the dimension of `m_coord`, `m_val` and
       * `m_normVal`.
       * */
      int m_numproj;

      /**
       * The i-th component of this vector contains a string storing the 
       * components used to get the i-th figure of merit of m_normVal.
       * */
      std::vector<std::string> m_coord;
  };
}

#endif

#include "latmrg/MeritProj.h"

#include "latticetester/Util.h"

#include <iomanip>
#include <sstream>

namespace LatMRG {

  MeritProj::MeritProj(int numproj) : Merit() {
    m_numproj = numproj;
    m_coord.reserve(numproj);
    setCoord();
    this->setDim(numproj);
  }

  //===========================================================================

  void MeritProj::setCoord() {
    for (unsigned int i = 0; i < m_coord.capacity(); i++) {
      m_coord.push_back(std::string());
    }
  }

  //===========================================================================

  std::string MeritProj::toString(bool rac, bool invert) {
    std::ostringstream os;
    double x[m_numproj], y[m_numproj];
    if (rac) {
      for (int i = 0; i < m_numproj; i++) {
        x[i] = LatticeTester::mysqrt(getMerit(i));
        y[i] = LatticeTester::mysqrt(getNormVal(i));
      }
    } else {
      for (int i = 0; i < m_numproj; i++) {
        x[i] = getMerit(i);
        y[i] = getNormVal(i);
      }
    }
    if (invert) {
      for (int i = 0; i < m_numproj; i++) {
        x[i] = 1.0/x[i];
      }
    }
    os << std::setprecision(5);
    os << "Projection:\t\t\tl_t\t\t\tS(I)\n";
    for (int i = 0; i < m_numproj; i++) {
      os << m_coord[i] << "\t" << x[i] << "\t" << y[i];
      os << std::endl;
    }
    return os.str();
  }
}

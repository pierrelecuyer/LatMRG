#include "latmrg/Projections.h"

namespace LatMRG {
  Projections::Projections(int dimProj, int minDim, std::vector<std::size_t>& projDim){
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

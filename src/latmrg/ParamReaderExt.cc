/**
  ParamReaderExt.cc for ISO C++
  version
  modified 28/05/06 10:36
authors: Hicham Wahbi
Frederik Rozon
Richard Simard
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <cstdlib>

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latmrg/ParamReaderExt.h"

using namespace std;


namespace LatMRG
{

  //===========================================================================

  ParamReaderExt::ParamReaderExt()
  {
    m_fileName.reserve(MAX_WORD_SIZE);
  }


  //===========================================================================

  ParamReaderExt::ParamReaderExt(std::string fileName) :ParamReader(fileName)
  {}



  //===========================================================================

  ParamReaderExt::~ParamReaderExt()
  {
    for (int i = 0; i < (int) m_lines.size(); i++)
      m_lines[i].clear();
    m_lines.clear();
  }


  //===========================================================================

  void ParamReaderExt::readGenType(GenType& field, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (strcasecmp(val.c_str(), "LCG") == 0)
      field = LCG;
    else if (strcasecmp(val.c_str(), "MRG") == 0)
      field = MRG;
    else if (strcasecmp(val.c_str(), "MWC") == 0)
      field = MWC;
    else if (strcasecmp(val.c_str(), "KOROBOV") == 0)
      field = KOROBOV;
    else if (strcasecmp(val.c_str(), "RANK1") == 0)
      field = RANK1;
    else if (strcasecmp(val.c_str(), "MMRG") == 0)
      field = MMRG;
    else
      LatticeTester::MyExit(1, "readGenType:   NO SUCH CASE");
  }


  //===========================================================================

  void ParamReaderExt::readLacunaryType(LacunaryType& lacunaryType, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (strcasecmp(val.c_str(), "NONE") == 0)
      lacunaryType = NONE;
    else if (strcasecmp(val.c_str(), "SUBVECTOR") == 0)
      lacunaryType = SUBVECTOR;
    else if (strcasecmp(val.c_str(), "ARBITRARYINDICES") == 0)
      lacunaryType = ARBITRARYINDICES;
    else
      LatticeTester::MyExit(1, "readLacunaryType:   NO SUCH CASE");
  }


  //===========================================================================

  void ParamReaderExt::readMMat(MMat & fields, unsigned int & ln, unsigned int pos,
      unsigned int numPos)
  {
    for (unsigned int i = pos; i < numPos; i++){
      for (unsigned int j = pos; j < numPos; j++){
        readMScal(fields(i,j), ln, j);
      }
      if(i!=numPos-1)
        ln++;
    }

  }


  //===========================================================================

  void ParamReaderExt::readInterval (MVect & B, MVect & C, unsigned int & ln, int k)
  {
    long m1, m2, m3;
    for (int i = 0; i < k; i++) {
      readNumber3 (B[i], m1, m2, m3, ++ln, 0);
      readNumber3 (C[i], m1, m2, m3, ++ln, 0);
      assert (C[i] >= B[i]);
    }
  }


  //===========================================================================

  void ParamReaderExt::readCalcType (LatticeTester::CalcType & field, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (0 == strcasecmp(val.c_str(), "PAL"))
      field = LatticeTester::PAL;
    else if (0 == strcasecmp(val.c_str(), "NORMPAL"))
      field = LatticeTester::NORMPAL;
    else if (0 == strcasecmp(val.c_str(), "BAL"))
      field = LatticeTester::BAL;
    else if (0 == strcasecmp(val.c_str(), "SEEKPAL"))
      field = LatticeTester::SEEKPAL;
    else
      LatticeTester::MyExit(1, "readCalcType:   NO SUCH CASE");
  }


  //===========================================================================

  void ParamReaderExt::readDecompType (DecompType & field, unsigned int line,
      unsigned int pos)
  {
    string val;
    getToken(val, line, pos);

    if (0 == strcasecmp(val.c_str(), "decomp"))
      field = DECOMP;
    else if (0 == strcasecmp(val.c_str(), "write"))
      field = DECOMP_WRITE;
    else if (0 == strcasecmp(val.c_str(), "read"))
      field = DECOMP_READ;
    else if (0 == strcasecmp(val.c_str(), "prime"))
      field = DECOMP_PRIME;
    else
      LatticeTester::MyExit(1, "readDecompType:   NO SUCH CASE");
  }


  //===========================================================================

  void ParamReaderExt::readLatticeType(LatticeType& field, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (0 == strcasecmp(val.c_str(), "FULL"))
      field = FULL;
    else if (0 == strcasecmp(val.c_str(), "RECURRENT"))
      field = RECURRENT;
    else if (0 == strcasecmp(val.c_str(), "ORBIT"))
      field = ORBIT;
    else if (0 == strcasecmp(val.c_str(), "PRIMEPOWER"))
      field = PRIMEPOWER;
    else
      LatticeTester::MyExit(1, "readLatticeType:   NO SUCH CASE");
  }


  //===========================================================================

  void ParamReaderExt::readOrbit (int J, MRGComponent<MScal> **comp, unsigned int & ln)
  {
    for (int j = 0; j < J; j++) {
      unsigned int k = comp[j]->k;
      readMVect (comp[j]->orbitSeed, ln, 1U, k, 1);
      ln++;
    }
  }


  //===========================================================================

  void ParamReaderExt::readImplemCond(ImplemCond& field, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (0 == strcasecmp(val.c_str(), "NoCond"))
      field = NO_COND;
    else if (0 == strcasecmp(val.c_str(), "AppFact"))
      field = APP_FACT;
    else if (0 == strcasecmp(val.c_str(), "NonZero")) {
      field = ZERO_COEF;
    } else if (0 == strcasecmp(val.c_str(), "Equal"))
      field = EQUAL_COEF;
    else if (0 == strcasecmp(val.c_str(), "PowerTwo"))
      field = POWER_TWO;
    else
      LatticeTester::MyExit(1, "readImplemCond:   NO SUCH CASE");
  }

  //===========================================================================

  void ParamReaderExt::readSearchMethod(SearchMethod& field, unsigned int ln, unsigned int pos)
  {
    string val;
    getToken(val, ln, pos);

    if (0 == strcasecmp(val.c_str(), "Exhaust"))
      field = EXHAUST;
    else if (0 == strcasecmp(val.c_str(), "Random"))
      field = RANDOM;
    else
      LatticeTester::MyExit(1, "readSearchMethod:   NO SUCH CASE");
  }

  //===========================================================================

  void ParamReaderExt::readLacunary(int ordre, int fromDim, int toDim,
      unsigned int & ln, bool & lacunary, int & lacGroupSize, BScal & lacSpacing,
      BVect & Lac, GenType genType)
  {
    readInt (lacGroupSize, ++ln, 0);
    const int t = lacGroupSize;
    if (t > 0)
      readZZ(lacSpacing, ln, 1);

    if (((t == 1) && (lacSpacing == 1)) || (t > toDim)) {
      lacunary = false;
      if (RANK1 == genType)
        return;
      if (toDim <= ordre)
        LatticeTester::MyExit(2, "ParamReaderExt::ReadLacunary:   toDim <= k");
      return;
    }

    lacunary = true;
    LatticeTester::CreateVect (Lac, toDim-1);

    int i;
    if (t < 0) {
      for (i = 0; i < toDim; i++)
        readBScal (Lac[i], ++ln, 0);
      //cout << "BVect Lac = " << Lac << endl;
      return;
    }

    NTL::ZZ Q1, Q;
    if (t > 0)
      Q1 = lacSpacing;

    Q = 1;
    // before Q was starting at 1, but this was producing a shift in the
    // indices of the lacunary vector.
    i = 0;
    while (true) {
      for (int j = 0; j < t; j++) {
        if (i < toDim) {
          conv (Lac[i], Q + j);
          i++;
        } else 
          return;
      }
      Q += Q1;
    }
  }

  //===========================================================================

  void ParamReaderExt::readMMRGLacunary(int ordre, int fromDim, int toDim,
      unsigned int & ln, bool & lacunary, LacunaryType & lacunaryType, int & numberLacIndices,
      BVect & Lac, GenType genType)
  {
    readLacunaryType(lacunaryType, ln, 0);

    if (lacunaryType == SUBVECTOR) {

      readInt (numberLacIndices, ln, 1);
      lacunary = true;

      // storing the lacunary indices for each vector
      BVect vectorSubLac;
      LatticeTester::CreateVect (vectorSubLac, numberLacIndices-1); 
      for (int i = 0; i < numberLacIndices; i++) {
        if (vectorSubLac[i] > ordre)
          LatticeTester::MyExit(1, "Lacunary indice too large. Must be smaller than the size of the generated vectors");
        readBScal (vectorSubLac[i], ++ln, 0);
      }

      // using those lacunary indices for each vector to build the arbitrary lacunary indices
      LatticeTester::CreateVect (Lac, toDim-1);
      int alpha = 0;
      for (int i = 0; i < toDim; i++) {
        Lac[i] = vectorSubLac[i - alpha * numberLacIndices] + alpha * ordre;
        if ( (i+1) % numberLacIndices == 0 )
          alpha++;
      } 
      return;

    } else if (lacunaryType == ARBITRARYINDICES) {

      readInt (numberLacIndices, ln, 1);
      lacunary = true;
      LatticeTester::CreateVect (Lac, numberLacIndices-1);
      for (int i = 0; i < numberLacIndices; i++)
        readBScal (Lac[i], ++ln, 0);
      return;

    } else if (lacunaryType == NONE) {

      lacunary = false;
      return;
    }
  }

  //===========================================================================

  bool ParamReaderExt::checkBound (const MScal & m, const MVect & A, int k)
  {
    for (int i = 0; i < k; i++) {
      assert (A[i] < m);
      assert (A[i] > -m);
    }
    return true;
  }

  //===========================================================================

  bool ParamReaderExt::checkBound (const MScal & m, const MMat & A, int k)
  {
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        assert (A[i][j] < m);
        assert (A[i][j] > -m);
      }
    }
    return true;
  }


  //===========================================================================

  bool ParamReaderExt::checkPrimePower(LatticeType lat, long m2, long m3, int k)
  {
    if (lat != PRIMEPOWER)
      return true;
    if (m2 <= 0)
      LatticeTester::MyExit(1, "LatticeType = PRIMEPOWER:   e <= 0");
    if (m3 != 0)
      LatticeTester::MyExit(1, "LatticeType = PRIMEPOWER:   c != 0");
    if (k > 1)
      LatticeTester::MyExit(1, "LatticeType = PRIMEPOWER:   k > 1");
    return true;
  }


  //===========================================================================

} // End namespace LatMRG

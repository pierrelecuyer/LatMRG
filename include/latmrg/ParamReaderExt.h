#ifndef LATMRG_PARAMREADEREXT_H
#define LATMRG_PARAMREADEREXT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <vector>


#include "NTL/ZZ.h"

#include "latticetester/Const.h"
#include "latticetester/ParamReader.h"
#include "latticetester/Util.h"

#include "latmrg/Const.h"
#include "latmrg/MRGComponent.h"



namespace LatMRG {

  /**
   * Utility class used to read basic parameter fields in a configuration file.
   * Lines whose first non-blank character is a <tt>#</tt> are considered as
   * comments and discarded.
   *
   */
  template <typename Int, typename Dbl>
    class ParamReaderExt: public LatticeTester::ParamReader<Int, Int, Dbl> {
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;
      public:
        static const int MAX_WORD_SIZE = 64;

        /**
         * Constructor.
         */
        ParamReaderExt();

        /**
         * Constructor. Opens the file `fileName`.
         */
        ParamReaderExt (std::string fileName);

        /**
         * Destructor.
         */
        ~ParamReaderExt();

        /**
         * Reads a `GenType` from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readGenType (GenType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a `IntMat` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
         * line into `field`.
         */
        void readMMat (IntMat & fields, unsigned int & ln, unsigned int pos,
            unsigned int numPos);

        /**
         * Reads `2k Int` tokens into vectors `B` and `C`, starting at the
         * <tt>ln</tt>-th line. These represent a box \f$[B_i, C_i]\f$, \f$i =
         * 1, 2, …, k\f$. The \f$B_i, C_i\f$ must be given in the order \f$B_1,
         * C_1, B_2, C_2, …, B_k, C_k\f$, each on a line of its own. Each
         * coefficient may be given in the form described in `readNumber3`
         * above.
         */
        void readInterval (IntVec & B, IntVec & C, unsigned int & ln, int k);

        /**
         * Reads the type of calculation to do (<tt>PAL</tt> or <tt>BAL</tt>)
         * for the `PALPHA` test from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readCalcType (LatticeTester::CalcType & field, unsigned int line,
            unsigned int pos);

        /**
         * Reads the decomposition type from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readDecompType (DecompType & field, unsigned int line,
            unsigned int pos);

        /**
         * Reads a lattice type from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readLatticeType (LatticeType & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads the initial state for an orbit of a combined generator with
         * \f$J\f$ components `comp`, starting at the <tt>ln</tt>-th line. Each
         * line must have a set of \f$k_j\f$ numbers which are the seeds of the
         * orbit for each component, where \f$k_j\f$ is the order of the
         * \f$j\f$-th component. On exiting this method, the line number is
         * reset to `ln` + \f$J\f$.
         */
        void readOrbit (int J, MRGComponent<Int> **comp, unsigned int & ln);

        /**
         * In the case where `lat = PRIMEPOWER`, checks that the modulus
         * \f$m\f$ is given in the form \f$m=b^e\f$, (that is \f$e>0\f$ and
         * \f$c=0\f$. Checks also that the order \f$k=1\f$. If these conditions
         * are not satisfied, stops the program.
         */
        bool checkPrimePower (LatticeType lat, long e, long c, int k);

        /**
         * Reads the fields, starting at the <tt>ln</tt>-th line, for a
         * generator of order \f$k\f$ and for a test in dimensions `fromDim` to
         * `toDim`, to determine the lacunary indices, which are positive
         * integers that will be read into `Lac`. The vector of lacunary
         * indices `Lac` is created here with the appropriate dimension. If the
         * first two values read are \f$s = \f$ `lacGroupSize` and \f$d = \f$
         * `lacSpacing` and if \f$s>0\f$, then what will be analyzed is the
         * lattice structure of vectors of the form \f$(u_{i+1}, …, u_{i+s},
         * u_{i+d+1},…, u_{i+d+s}, u_{i+2d+1},…, u_{i+2d+s}, …)\f$, formed by
         * groups of \f$s\f$ successive values, taken \f$d\f$ values apart. To
         * analyze lacunary indices that are not evenly spaced, put \f$s
         * = -t\f$ where <em>\f$t = \f$ MaxDim</em> and, on the \f$t\f$ lines
         * that follow, give the \f$t\f$ lacunary indices \f$i_1,…,i_t\f$,
         * which are to be interpreted as in Section. In all these cases, the `lacunary`
         * flag is set `true`. To analyze vectors of successive values (the
         * non-lacunary case), take \f$s=d=1\f$ or \f$s\ge\f$ *MaxDim*. In
         * this case, the `lacunary` flag is set `false`.
         */
        void readLacunary (int k, int fromDim, int toDim, unsigned int & ln,
            bool & lacunary, int & lacGroupSize, Int & lacSpacing,
            IntVec & Lac, GenType genType);

        /**
         * Reads the lacunary configuration specified in the input file. If
         * *lacunaryType* is NONE, then the projection is done without lacunary
         * indices. If *lacunaryType* is SUBVECTOR, the projection is performed
         * using only the specified coordinates for each vector (*m_numberLacIndices*
         * coordinates per vector) until the dimension *toDim* is reached. Finally,
         * if *lacunaryType* is ARBITRARYINDICES, the projection is performed using
         * the specified indices (*m_numberLacIndices* of them).
         */
        void readMMRGLacunary(int ordre, int fromDim, int toDim, unsigned int & ln,
            bool & lacunary, LacunaryType & lacunaryType, int & m_numberLacIndices,
            IntVec & Lac, GenType genType);

        /**
         * Reads the type of lacunary projection used (either NONE, SUBVECTOR, or
         * ARBITRARYINDICES).
         */
        void readLacunaryType(LacunaryType& lacunaryType, unsigned int ln, unsigned int pos);

        /**
         * Reads an implementation condition from the <tt>pos</tt>-th token of
         * the <tt>ln</tt>-th line into `field`.
         */
        void readImplemCond (ImplemCond & field, unsigned int ln,
            unsigned int pos);

        /**
         * Reads a search method from the <tt>pos</tt>-th token of the
         * <tt>ln</tt>-th line into `field`.
         */
        void readSearchMethod (SearchMethod & field, unsigned int ln,
            unsigned int pos);

        /**
         * Checks that the components of \f$A\f$ satisfy \f$-m < A_i < m\f$,
         * for \f$i=1, 2, …, k\f$.
         */
        bool checkBound (const Int & m, const IntVec & A, int k);

        /**
         * Checks that the components of \f$A\f$ satisfy \f$-m < A_{i,j} < m\f$,
         * for \f$i,j=1, 2, …, k\f$.
         */
        bool checkBound (const Int & m, const IntMat & A, int k);

      private:

        /**
         * Internal line buffer.
         */
        std::vector<std::string> m_lines;

        /**
         * The path of the opened file.
         */
        std::string m_fileName;

        /**
         * Checks if the character `c` is to be considered as a token separator
         * or not.
         */
        bool isDelim (char c);

        /**
         * Does nothing for now.
         */
        void init() {}
    };

  template <typename Int, typename Dbl>
    ParamReaderExt<Int, Dbl>::ParamReaderExt()
    {
      m_fileName.reserve(MAX_WORD_SIZE);
    }


  //===========================================================================

  template <typename Int, typename Dbl>
    ParamReaderExt<Int, Dbl>::ParamReaderExt(std::string fileName):
      LatticeTester::ParamReader<Int, Int, Dbl>(fileName)
  {}



  //===========================================================================

  template <typename Int, typename Dbl>
    ParamReaderExt<Int, Dbl>::~ParamReaderExt()
    {
      for (int i = 0; i < (int) m_lines.size(); i++)
        m_lines[i].clear();
      m_lines.clear();
    }

  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readGenType(GenType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

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
      else if (strcasecmp(val.c_str(), "COMBO") == 0)
        field = COMBO;
      else
        LatticeTester::MyExit(1, "readGenType:   NO SUCH CASE");
    }

  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readLacunaryType(LacunaryType& lacunaryType, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readMMat(IntMat & fields, unsigned int & ln, unsigned int pos,
        unsigned int numPos)
    {
      for (unsigned int i = pos; i < numPos; i++){
        for (unsigned int j = pos; j < numPos; j++){
          this->readMScal(fields(i,j), ln, j);
        }
        if(i!=numPos-1)
          ln++;
      }

    }


  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readInterval (IntVec & B, IntVec & C, unsigned int & ln, int k)
    {
      long m1, m2, m3;
      for (int i = 0; i < k; i++) {
        this->readNumber3 (B[i], m1, m2, m3, ++ln, 0);
        this->readNumber3 (C[i], m1, m2, m3, ++ln, 0);
        assert (C[i] >= B[i]);
      }
    }


  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readCalcType (LatticeTester::CalcType & field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readDecompType (DecompType & field, unsigned int line,
        unsigned int pos)
    {
      std::string val;
      this->getToken(val, line, pos);

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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readLatticeType(LatticeType& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readOrbit (int J, MRGComponent<Int> **comp, unsigned int & ln)
    {
      for (int j = 0; j < J; j++) {
        unsigned int k = comp[j]->getK();
        this->readMVect (comp[j]->getOrbitSeed(), ln, 1U, k, 1);
        ln++;
      }
    }


  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readImplemCond(ImplemCond& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readSearchMethod(SearchMethod& field, unsigned int ln, unsigned int pos)
    {
      std::string val;
      this->getToken(val, ln, pos);

      if (0 == strcasecmp(val.c_str(), "Exhaust"))
        field = EXHAUST;
      else if (0 == strcasecmp(val.c_str(), "Random"))
        field = RANDOM;
      else
        LatticeTester::MyExit(1, "readSearchMethod:   NO SUCH CASE");
    }

  //===========================================================================

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readLacunary(int ordre, int fromDim, int toDim,
        unsigned int & ln, bool & lacunary, int & lacGroupSize, Int & lacSpacing,
        IntVec & Lac, GenType genType)
    {
      this->readInt (lacGroupSize, ++ln, 0);
      const int t = lacGroupSize;
      if (t > 0)
        this->readBScal(lacSpacing, ln, 1);

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
          this->readBScal (Lac[i], ++ln, 0);
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

  template <typename Int, typename Dbl>
    void ParamReaderExt<Int, Dbl>::readMMRGLacunary(int ordre, int fromDim, int toDim,
        unsigned int & ln, bool & lacunary, LacunaryType & lacunaryType, int & numberLacIndices,
        IntVec & Lac, GenType genType)
    {
      readLacunaryType(lacunaryType, ln, 0);

      if (lacunaryType == SUBVECTOR) {

        this->readInt (numberLacIndices, ln, 1);
        lacunary = true;

        // storing the lacunary indices for each vector
        IntVec vectorSubLac;
        LatticeTester::CreateVect (vectorSubLac, numberLacIndices-1);
        for (int i = 0; i < numberLacIndices; i++) {
          if (vectorSubLac[i] > ordre)
            LatticeTester::MyExit(1, "Lacunary indice too large. Must be smaller than the size of the generated vectors");
          this->readBScal (vectorSubLac[i], ++ln, 0);
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

        this->readInt (numberLacIndices, ln, 1);
        lacunary = true;
        LatticeTester::CreateVect (Lac, numberLacIndices-1);
        for (int i = 0; i < numberLacIndices; i++)
          this->readBScal (Lac[i], ++ln, 0);
        return;

      } else if (lacunaryType == NONE) {

        lacunary = false;
        return;
      }
    }

  //===========================================================================

  template <typename Int, typename Dbl>
    bool ParamReaderExt<Int, Dbl>::checkBound (const Int & m, const IntVec & A, int k)
    {
      for (int i = 0; i < k; i++) {
        assert (A[i] < m);
        assert (A[i] > -m);
      }
      return true;
    }

  //===========================================================================

  template <typename Int, typename Dbl>
    bool ParamReaderExt<Int, Dbl>::checkBound (const Int & m, const IntMat & A, int k)
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

  template <typename Int, typename Dbl>
    bool ParamReaderExt<Int, Dbl>::checkPrimePower(LatticeType lat, long m2, long m3, int k)
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

  extern template class ParamReaderExt<std::int64_t, double>;
  extern template class ParamReaderExt<NTL::ZZ, double>;
  extern template class ParamReaderExt<NTL::ZZ, NTL::RR>;

}
#endif

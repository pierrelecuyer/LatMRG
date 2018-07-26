#ifndef PARAMREADEREXT_H
#define PARAMREADEREXT_H

#include "NTL/ZZ.h"

#include "latticetester/Types.h"
#include "latticetester/Const.h"
#include "latticetester/ParamReader.h"

#include "latmrg/Const.h"
#include "latmrg/MRGComponent.h"

#include <string>
#include <vector>


namespace LatMRG {

  /**
   * Utility class used to read basic parameter fields in a configuration file.
   * Lines whose first non-blank character is a <tt>#</tt> are considered as
   * comments and discarded.
   *
   */
  class ParamReaderExt: public LatticeTester::ParamReader<MScal, BScal, RScal> {
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
       * Reads a `BMat` from the <tt>pos</tt>-th token of the <tt>ln</tt>-th
       * line into `field`.
       */
      void readMMat (MMat & fields, unsigned int & ln, unsigned int pos,
          unsigned int numPos);

      /**
       * Reads `2k MScal` tokens into vectors `B` and `C`, starting at the
       * <tt>ln</tt>-th line. These represent a box \f$[B_i, C_i]\f$, \f$i =
       * 1, 2, …, k\f$. The \f$B_i, C_i\f$ must be given in the order \f$B_1,
       * C_1, B_2, C_2, …, B_k, C_k\f$, each on a line of its own. Each
       * coefficient may be given in the form described in `readNumber3`
       * above.
       */
      void readInterval (MVect & B, MVect & C, unsigned int & ln, int k);

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
      void readOrbit (int J, MRGComponent<MScal> **comp, unsigned int & ln);

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
       * which are to be interpreted as in Section  {@link
       * REF__sec1_sec_lacunary lacunary}. In all these cases, the `lacunary`
       * flag is set `true`. To analyze vectors of successive values (the
       * non-lacunary case), take \f$s=d=1\f$ or \f$s\ge\f$ *MaxDim*. In
       * this case, the `lacunary` flag is set `false`.
       */
      void readLacunary (int k, int fromDim, int toDim, unsigned int & ln,
          bool & lacunary, int & lacGroupSize, BScal & lacSpacing,
          BVect & Lac, GenType genType);

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
          BVect & Lac, GenType genType);

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
      bool checkBound (const MScal & m, const MVect & A, int k);

      /**
       * Checks that the components of \f$A\f$ satisfy \f$-m < A_{i,j} < m\f$,
       * for \f$i,j=1, 2, …, k\f$.
       */
      bool checkBound (const MScal & m, const MMat & A, int k);

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

}
#endif

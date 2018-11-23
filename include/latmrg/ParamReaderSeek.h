#ifndef PARAMREADERSEEK_H
#define PARAMREADERSEEK_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <string>

#include "latticetester/Util.h"

#include "latmrg/ParamReaderExt.h"
#include "latmrg/SeekConfig.h"

namespace LatMRG {

  /**
   * This class is used to read a configuration file for the executable
   * programs `seek*`, created from the `SeekMain` main program. The format of
   * the configuration file is described in this guide for the program
   * `SeekMain` on page (FIXME: page#).
   *
   */
  template<typename Int, typename Dbl>
    class ParamReaderSeek : public ParamReaderExt<Int, Dbl> {
      public:

        /**
         * Default Constructor.
         */
        ParamReaderSeek();

        /**
         * Constructor. `fname` is the name of the configuration file.
         */
        ParamReaderSeek (std::string fname);

        /**
         * Destructor.
         */
        ~ParamReaderSeek();

        /**
         * Variables initialization.
         */
        void init();

        /**
         * Reads the configuration parameters and puts them in an instance of
         * <tt>SeekConfig</tt>.
         */
        void read (SeekConfig<Int> & config);
    };

  //===========================================================================

  template<typename Int, typename Dbl>
    ParamReaderSeek<Int, Dbl>::ParamReaderSeek ()
    {}


  //===========================================================================

  template<typename Int, typename Dbl>
    ParamReaderSeek<Int, Dbl>::ParamReaderSeek (std::string fileName): ParamReaderExt<Int, Dbl> (fileName)
  {}


  //===========================================================================

  template<typename Int, typename Dbl>
    ParamReaderSeek<Int, Dbl>::~ParamReaderSeek ()
    {
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void ParamReaderSeek<Int, Dbl>::read (SeekConfig<Int> & config)
    {
      this->getLines ();
      int ln = 0;
      unsigned int lnu;
      int s;

      this->readCriterionType (config.criter, ln, 0);
      if (config.criter == LatticeTester::SPECTRAL)
        this->readNormaType (config.normaType, ln, 1);

      this->readBool (config.readGenFile, ++ln, 0);
      if (config.readGenFile) {
        this->readString (config.fileName, ln, 1);
        config.fileName += ".gen";
      }
      this->readInt (config.J, ++ln, 0);

      config.Ms = new Int[config.J];
      config.Ks = new int[config.J];

      config.compon.reserve (config.J);
      int maxOrder = 0;
      Component<Int> comp;
      for (int i = 0; i < config.J; i++) {
        this->readGenType (comp.genType, ++ln, 0);
        if (comp.genType == MMRG){
          LatticeTester::MyExit(1, "Not Yet Implemented");
        }
        this->readNumber3 (comp.modulus.m, comp.modulus.b, comp.modulus.e,
            comp.modulus.c, ++ln, 0);
        if (0 == comp.modulus.e)
          comp.modulus.init (comp.modulus.m);
        else
          comp.modulus.init (comp.modulus.b, comp.modulus.e, comp.modulus.c);
        config.Ms[i] = comp.modulus.m;

        this->readInt (comp.k, ++ln, 0);
        if (maxOrder < comp.k)
          maxOrder = comp.k;
        config.Ks[i] = comp.k;
        this->readBool (comp.PerMax, ++ln, 0);

        this->readDecompType (comp.F1, ++ln, 0);
        if (comp.F1 == DECOMP_WRITE || comp.F1 == DECOMP_READ) {
          comp.file1.reserve (this->MAX_WORD_SIZE);
          this->readString (comp.file1, ln, 1);
        }
        this->readDecompType (comp.F2, ++ln, 0);
        if (comp.F2 == DECOMP_WRITE || comp.F2 == DECOMP_READ) {
          comp.file2.reserve (this->MAX_WORD_SIZE);
          this->readString (comp.file2, ln, 1);
        }

        this->readImplemCond (comp.implemCond, ++ln, 0);
        if (POWER_TWO == comp.implemCond) {
          this->readInt (comp.NumBits, ln, 1);
          this->readInt (comp.HighestBit, ln, 2);
        } else if (EQUAL_COEF == comp.implemCond ||
            ZERO_COEF == comp.implemCond) {
          this->readInt (comp.ncoef, ln, 1);
          comp.Icoef = new int[1 + comp.ncoef];
          for (s = 1; s < comp.ncoef; s++)
            this->readInt (comp.Icoef[s], ln, 1 + s);
          comp.Icoef[0] = 0;
          comp.Icoef[comp.ncoef] = comp.k;
        }

        comp.b.SetLength (comp.k);
        comp.c.SetLength (comp.k);

        lnu = ln;
        if (ZERO_COEF == comp.implemCond) {
          LatticeTester::SetZero (comp.b, comp.k);
          LatticeTester::SetZero (comp.c, comp.k);
          this->readInterval (comp.b, comp.c, lnu, comp.ncoef);
          for (s = comp.ncoef-1; s >= 0; s--) {
            comp.b[comp.Icoef[s]] = comp.b[s];
            comp.c[comp.Icoef[s]] = comp.c[s];
            if (s < comp.Icoef[s]) {
              comp.b[s] = 0;
              comp.c[s] = 0;
            }
          }

        } else if (EQUAL_COEF == comp.implemCond) {
          this->readInterval (comp.b, comp.c, lnu, comp.ncoef);
          for (s = comp.ncoef-1; s >= 0; s--) {
            for (int i = comp.Icoef[s]; i > comp.Icoef[s-1]; i--) {
              comp.b[i] = comp.b[s];
              comp.c[i] = comp.c[s];
            }
          }
        } else
          this->readInterval (comp.b, comp.c, lnu, comp.k);
        ln = lnu;

        if (comp.genType == RANK1)
          comp.b[0] = comp.c[0] = 1;
        this->checkBound (comp.modulus.m, comp.b, comp.k);
        this->checkBound (comp.modulus.m, comp.c, comp.k);

        if (POWER_TWO == comp.implemCond) {
          for (int i = 0; i < comp.k; i++) {
            comp.b[i] = 0;
            comp.c[i] = comp.HighestBit;
          }
        }
        this->readSearchMethod (comp.searchMethod, ++ln, 0);
        if (RANDOM == comp.searchMethod) {
          this->readInt (comp.numReg, ln, 1);

          this->readInt (comp.H, ln, 2);
          this->readInt (comp.Hk, ln, 3);
        } else
          comp.numReg = 1;

        config.compon.push_back (comp);
      }

      this->readInt (config.C, ++ln, 0);
      this->readDoubleVect (config.minMerit, ++ln, 0, config.C, 0);
      this->readDoubleVect (config.maxMerit, ++ln, 0, config.C, 0);
      this->readIntVect (config.numGen, ++ln, 0, config.C, 0);
      this->readInt (config.d, ++ln, 0);
      if (config.d < 1)
        LatticeTester::MyExit (1, "ParamReaderSeek:   config.d < 1");
      config.td = new int[config.d + 1];
      this->readIntVect (config.td, ++ln, 0, config.d + 1, 0);
      if (config.td[0] <= maxOrder)
        config.td[0] = maxOrder + 1;

      this->readBool (config.dualF, ++ln, 0);
      this->readInt (config.speed, ++ln, 0);
      this->readLatticeType (config.latType, ++ln, 0);

      this->readInt (config.lacGroupSize, ++ln, 0);
      this->readInt (config.lacSpacing, ln, 1);
      this->readLong (config.maxNodesBB, ++ln, 0);

      double tem;
      this->readDouble (tem, ++ln, 0);
      char c;
      this->readChar (c, ln, 1);
      switch (c) {
        case 's':
          config.duration = tem;
          break;
        case 'm':
          config.duration = tem * 60;
          break;
        case 'h':
          config.duration = tem * 3600;
          break;
        case 'd':
          config.duration = tem * 86400;
          break;
        default:
          config.duration = 0;
          std::cout << "***** reading duration:  IMPOSSIBLE CASE";
      }

      this->readLong (config.seed, ++ln, 0); // seed of the generator
      //   this->readLong (config.s2, ln, 2);
      this->readOutputType (config.outputType, ++ln, 0);
    }

}
#endif

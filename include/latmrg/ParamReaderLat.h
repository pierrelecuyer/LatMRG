#ifndef PARAMREADERLAT_H
#define PARAMREADERLAT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstring>

#include "latticetester/Util.h"
#include "latticetester/Const.h"

#include "latmrg/Const.h"
#include "latmrg/ParamReaderLat.h"
#include "latmrg/MRGComponentFactory.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/LatConfig.h"


namespace LatMRG {

  /**
   * This class is used to read a configuration file for the executable
   * programs `lat*`, created from the `LatMain` main program. The format of
   * the configuration file is described in this guide for the program
   * `LatMain` on page (FIXME: page#).
   *
   */
  template<typename Int, typename Dbl>
    class ParamReaderLat : public ParamReaderExt<Int, Dbl> {
      public:

        /**
         * Constructor. Opens configuration file `fileName`.
         */
        ParamReaderLat (std::string fileName): ParamReaderExt<Int, Dbl> (fileName)
      {}

        /**
         * Destructor.
         */
        ~ParamReaderLat() {}

        /**
         * Reads the configuration file into `config` for the Beyer and the
         * spectral tests.
         */
        void read (LatConfig<Int> & config)
        {
          this->getLines ();
          unsigned int ln = 0;

          this->readCriterionType (config.criter, ln, 0);
          if (config.criter == LatticeTester::SPECTRAL)
            this->readNormaType (config.norma, ln, 1);
          if (config.criter == LatticeTester::PALPHA) {
            this->readCalcType (config.calcPalpha, ++ln, 0);
            config.J = 1;
          } else {
            this->readNormType (config.norm, ++ln, 0);
            this->readBool (config.readGenFile, ++ln, 0);
            if (config.readGenFile)
              this->readString (config.fileName, ln, 1);
            this->readInt (config.J, ++ln, 0); // J
          }

          config.setJ (config.J);

          Int m;
          long m1(0), m2(0), m3(0);
          int k(0);
          NTL::vector<Int> a;
          NTL::matrix<Int> A;
          int maxOrder = 0;
          std::string coefkind (" ");  // = NOCOND, EQUAL, NONZERO
          for (int i = 0; i < config.J; i++) {
            this->readGenType (config.genType[i], ++ln, 0);
            this->readNumber3 (m, m1, m2, m3, ++ln, 0);
            if (config.criter == LatticeTester::PALPHA)
              k = 1;
            else
              this->readInt (k, ++ln, 0);
            if (maxOrder < k)
              maxOrder = k;
            if (config.genType[i] == MMRG){

              A.resize(k, k);
              this->readMMat(A, ++ln, 0, k);
              this->checkBound(m, A, k);
            }
            else{
              a.SetLength (k);
              LatticeTester::SetZero (a, k);

              if (config.criter != LatticeTester::PALPHA) {
                this->readString (coefkind, ++ln, 0);
              }
              if (0 == strcasecmp(coefkind.c_str(), "NONZERO")) {
                int s;
                this->readInt (s, ln, 1);
                int* T = new int[s+3];
                T[0] = 0;
                this->readIntVect (T, ln++, 2, s+3, 0);
                this->readMVect (a, ln, 1, s, 0);
                //   ln++;
                for (int j = s; j >= 1; j--) {
                  int r = T[j];
                  a[r] = a[j-1];
                  a[j-1] = 0;
                }
                delete[] T;
              } else if (0 == strcasecmp(coefkind.c_str(), "EQUAL")) {
                int s;
                this->readInt (s, ln, 2);
                int* T = new int[s+3];
                this->readIntVect (T, ln++, 3, s+3, 1);
                int p0 = 1;
                for (int j = 1; j <= s; j++) {
                  int r = T[j];
                  this->readMScal (a[r], ln, j);
                  for (int p = p0; p < r; p++) {
                    a[p] = a[r];
                  }
                  p0 = r + 1;
                }
                delete[] T;
              } else   // NoCond
                this->readMVect (a, ++ln, 0, k, 0);
              this->checkBound (m, a, k);
            }

            if (config.genType[i] == KOROBOV) {
              if (1 != k)
                LatticeTester::MyExit(1, "KOROBOV must have k = 1");
              config.comp[i] = new MRGComponent<Int> (m, a, 1);
            } else if (config.genType[i] == RANK1) {
              config.comp[i] = new MRGComponent<Int> (m, a, k);
            } else if (config.genType[i] == MRG || config.genType[i] == LCG) {
              config.comp[i] = new MRGComponent<Int> (m1, m2, m3, a, k);
              if (1 == k)
                config.comp[i]->module.reduceM (config.comp[i]->a[0]);
            } else if (config.genType[i] == MWC) {
              config.comp[i] = MRGComponentFactory<Int>::fromMWC (m, a, k);
            }
            else if (config.genType[i] == MMRG){
              config.comp[i] = new MRGComponent<Int>(m, A, k);
            } else
              assert(0);
          }

          this->readInt (config.d, ++ln, 0);
          if (config.d < 1)
            LatticeTester::MyExit (1, "ParamReaderLat:   config.d < 1");
          config.td = new int[config.d+1];
          this->readIntVect (config.td, ++ln, 0, 1 + config.d, 0);
          int fromDim = config.td[0];
          int toDim = config.td[1];

          if (config.criter == LatticeTester::PALPHA) {
            this->readBool (config.primeM, ++ln, 0);
            this->readBool (config.verifyM, ln, 1);
            this->readBool (config.maxPeriod, ++ln, 0);
            this->readBool (config.verifyP, ln, 1);
            this->readInt (config.alpha, ++ln, 0);
            this->readInt (config.seed, ++ln, 0);
            config.Beta.reserve(1 + toDim);
            ++ln;
            for (int i = 0; i <= toDim; i++)
              this->readDouble (config.Beta[i], ln, 1 + i);
            config.lacunary = false;

          } else {
            this->readBool (config.dualF, ++ln, 0);
            this->readLatticeType (config.latType, ++ln, 0);
            if ((config.genType[0] == RANK1 || config.genType[0] == KOROBOV) &&
                config.latType != FULL) {
              LatticeTester::MyExit (1,
                  "ParamReaderLat:   latType must be FULL for KOROBOV or RANK1 lattices");
            }
            this->checkPrimePower (config.latType, m2, m3, k);
            if (config.latType == ORBIT) {
              LatticeTester::MyExit (1, "case ORBIT is not finished");
              this->readOrbit (config.J, config.comp, ++ln);
              --ln;
            }

            if (config.genType[0] == MMRG) {
              this->readMMRGLacunary (k, fromDim, toDim, ++ln, config.lacunary, config.lacunaryType, 
                  config.numberLacIndices, config.Lac, config.genType[0]);
            } else {
              this->readLacunary (k, fromDim, toDim, ln, config.lacunary, config.lacGroupSize, 
                  config.lacSpacing, config.Lac, config.genType[0]);
            }


            if (config.genType[0] != RANK1 && config.td[0] <= maxOrder &&
                !config.lacunary)
              config.td[0] = maxOrder + 1;

            this->readLong (config.maxNodesBB, ++ln, 0);
            this->readBool (config.invertF, ++ln, 0);
            this->readInt (config.detailF, ++ln, 0);
          }
          this->readOutputType (config.outputType, ++ln, 0);
        }
    };

}
#endif

#ifndef LATMRGROUTINES_H
#define LATMRGROUTINES_H

#include "latticetester/Rank1Lattice.h"
#include "latticetester/NormaPalpha.h"

#include "latmrg/MMRGLattice.h"
#include "latmrg/KorobovLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/TestProjections.h"

namespace LatMRG {
  /** \file LatMRGRoutines.h
   * This file contains high level routines that allow the user to performed all
   * the calculations specified by a LatConfig object.
   * */

  /**
   * This function will compute the figure of merit for all the projections
   * required by `config`. `config` will contains a lot of information that we
   * will not cover in depth here, see *LatConfig* for more info.
   */
  template<typename Int, typename BasInt, typename BasVec, typename BasMat,
    typename Dbl, typename DblVec, typename RedDbl>
      std::vector<double> ComputeFigureOfMerit (LatConfig<Int>& config)
      {

        //Writer* rw = createWriter (infile, config.outputType);

        std::vector<double> result;

        LatticeTester::IntLattice<Int, BasInt, Dbl, RedDbl> *lattice = 0;
        LatticeTester::IntLattice<Int, BasInt, Dbl, RedDbl> *master = 0;
        LatticeTester::Lacunary<BasInt> *plac = 0;
        int toDim = config.td[1];
        int fromDim = config.td[0];
        bool memLacF = true; // Lacunary with only used lines-columns of bases
        //memLacF = false; // Lacunary with all lines-columns of bases

        if (config.J > 1) { //Several MRG
          lattice = MRGLatticeFactory<Int>::fromCombMRG (config.comp, config.J,
              toDim, 0, config.latType, config.norm);

        } else {
          if (config.latType == PRIMEPOWER) {
            config.comp[0]->module.reduceM (config.comp[0]->a[0]);
            if (memLacF && config.lacunary)
              lattice = new MRGLatticeLac<Int> (config.comp[0]->module.mRed,
                  config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                  config.latType, config.norm);
            else{
              lattice = new MRGLattice<Int> (config.comp[0]->module.mRed,
                  config.comp[0]->a, toDim, config.comp[0]->k,
                  config.latType, config.norm);

            }

          } else if (config.genType[0] == MRG || config.genType[0] == LCG) {

            if (memLacF && config.lacunary){            
              lattice = new MRGLatticeLac<Int> (config.comp[0]->module.mRed,
                  config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                  config.latType, config.norm);
            }
            else{
              lattice = new MRGLattice<Int> (config.comp[0]->module.mRed,
                  config.comp[0]->a, toDim, config.comp[0]->k,
                  config.latType, config.norm);
            }

          } else if (config.genType[0] == KOROBOV) {
            lattice = new KorobovLattice<Int> (config.comp[0]->getM (),
                config.comp[0]->a[1], toDim, config.norm);
          } else if (config.genType[0] == RANK1) {
            lattice = new LatticeTester::Rank1Lattice<Int, BasInt, Dbl, RedDbl>
              (config.comp[0]->getM (), config.comp[0]->a, config.comp[0]->k,
               config.norm);

          } else if (config.genType[0] == MMRG) {
            if (memLacF && config.lacunary) {
              lattice = new MMRGLattice<Int> (config.comp[0]->getM(),
                  config.comp[0]->A, toDim,config.comp[0]->k, config.lacunaryType,
                  config.Lac, config.norm);           
            } else {
              lattice = new MMRGLattice<Int> (config.comp[0]->getM(), 
                  config.comp[0]->A, toDim,config.comp[0]->k, config.norm);
            }
          }
        }

        double minVal[1 + toDim];
        LatticeTester::SetZero (minVal, toDim);

        LatticeTester::Normalizer<RedDbl> *normal = 0;

        if (config.criter == LatticeTester::SPECTRAL) {
          normal = lattice->getNormalizer (config.norma, 0, config.dualF);
          // creates and returns the normalizer corresponding to config.norma
          normal->setNorm (config.norm);
        } else if (config.criter == LatticeTester::PALPHA &&
            (config.calcPalpha == LatticeTester::NORMPAL ||
             config.calcPalpha == LatticeTester::BAL)) {
          normal = new LatticeTester::NormaPalpha<Int, RedDbl> (
              lattice->getModulo(), config.alpha, toDim);
        }

        if (!memLacF && config.lacunary) {
          plac = new LatticeTester::Lacunary<BasInt> (config.Lac, toDim);
          lattice->setLac (*plac);
        }

        switch (config.criter) {
          case LatticeTester::SPECTRAL: 
            {
              LatTestSpectral spectralTest (normal, lattice);
              lattice->buildBasis (fromDim - 1);
              //spectralTest.attach (&report);
              //report.printHeader ();
              spectralTest.setDualFlag (config.dualF);
              spectralTest.setInvertFlag (config.invertF);
              spectralTest.setDetailFlag (config.detailF);
              spectralTest.setMaxAllDimFlag (true);
              spectralTest.setMaxNodesBB (config.maxNodesBB);

              if (1 == config.d) {
                spectralTest.test (fromDim, toDim, minVal);
                lattice->write();
                //footer.setLatticeTest (&spectralTest);
                //report.printTable ();
                //report.printFooter ();
              } else {
                if (config.genType[0] == MRG || config.genType[0] == LCG)
                  master = new MRGLattice<Int> (*(MRGLattice<Int> *) lattice);
                else if (config.genType[0] == KOROBOV)
                  master = new KorobovLattice<Int> (
                      *(KorobovLattice<Int> *) lattice);
                else if (config.genType[0] == RANK1)
                  master = new LatticeTester::Rank1Lattice<Int, BasInt, Dbl, RedDbl>
                    (*(LatticeTester::Rank1Lattice<Int, BasInt, Dbl, RedDbl> *) lattice);
                else if (config.genType[0] == MMRG)
                  master = new MMRGLattice<Int> (*(MMRGLattice<Int> *) lattice);
                else {
                  std::cout << "GenType has not been specified correctly in config.\n";
                }

                master->buildBasis (toDim);
                TestProjections proj (master, lattice, &spectralTest, config.td,
                    config.d);
                //proj. setOutput (rw);
                proj.setDualFlag (config.dualF);
                proj.setPrintF (false);
                double minVal[config.td[0]];
                if (config.genType[0] == MMRG)
                  proj.run(false, false, minVal);
                else
                  proj.run(true, false, minVal);
                /*rw->writeString ("\nMin merit:   ");
                  rw->writeDouble (sqrt (merit));
                  rw->newLine ();
                  rw->writeString ("Num projections:   ");
                  rw->writeInt (nbProj);
                  rw->newLine ();*/
                // nbProj = proj.calcNumProjections(stationary, false);
                // cout << "Num projections2:  " << nbProj << endl << endl;
                delete master;
              }

              // storing the figures of merit in result
              for (int i = fromDim; i <= toDim; i++) 
                result.push_back(spectralTest.getMerit().getNormVal(i));

            }
            break;

          case LatticeTester::BEYER: 
            {
              LatTestBeyer<Int, BasInt, BasVec, BasMat, Dbl, DblVec, RedDbl>
                beyerTest (lattice);
              lattice->buildBasis (fromDim - 1);
              //beyerTest.attach (&report);
              //report.printHeader ();
              beyerTest.setDualFlag (config.dualF);
              beyerTest.setMaxAllDimFlag (true);
              beyerTest.setMaxNodesBB (config.maxNodesBB);
              beyerTest.setDetailFlag (config.detailF);
              beyerTest.test (fromDim, toDim, minVal);
              //footer.setLatticeTest (&beyerTest);
              //report.printTable ();
              //report.printFooter ();
              //rw->writeString (lattice->toStringDualBasis ());

              // storing the figures of merit in result
              for (int i = fromDim; i <= toDim; i++) 
                result.push_back(beyerTest.getMerit().getNormVal(i));
            }
            break;

          case LatticeTester::PALPHA: 
            {
              LatTestPalpha<Int> palphaTest (normal, lattice);
              palphaTest.setConfig (&config);
              //palphaTest.attach (&report);
              //report.printHeader ();
              if (1 == config.d) {
                palphaTest.test (fromDim, toDim, minVal);
                //footer.setLatticeTest (&palphaTest);
                //report.printTable ();
                //report.printFooter ();
              } else {
                MRGLattice<Int> master = MRGLattice<Int> (*(MRGLattice<Int> *) 
                    lattice);
                master.buildBasis (toDim);
                TestProjections proj (&master, lattice, &palphaTest, config.td,
                    config.d);
                //proj. setOutput (rw);
                /*rw->writeString ("\n\nMin merit:   ");
                  rw->writeDouble (sqrt (merit));
                  rw->newLine ();
                  rw->writeString ("Num projections:   ");
                  rw->writeInt (nbProj);
                  rw->newLine ();
                  rw->newLine ();*/
              }

              // storing the figures of merit in result
              for (int i = fromDim; i <= toDim; i++) 
                result.push_back(palphaTest.getMerit().getNormVal(i));
            }
            break;

          default:
            std::cerr << "Default case for config.criter" << std::endl;
            return result;
        }

        if (normal != 0)
          delete normal;
        if (!memLacF && config.lacunary)
          delete plac;
        delete lattice;
        //delete rw;

        return result;
      }

  /**
   * Prints the figures of merit that have been computed and that are stored in 
   * result. This could be modified to just call the already existing functions
   * in the low level API from a LatConfig object instead.
   * */
  void printResult(const std::vector<double> & result, const int & fromDim);
  // ajouter plus de détails

  // une fonction qui calcule le shortest vector

  // une fonction d'aide à la configuration de LatConfig
  /** I have no idea what this function is supposed to do. It currently does 
   * nothing. I think the idea is to create a polymorphic function that 
   * automatically builds a LatConfig object ready to be tested with
   * ComputeFigureOfMerit.
   * */
  template<typename Int>
    void initConfigSpectralTest(LatConfig<Int>& config){}
  /*
   * parameters:
   * fromDIm, toDim
   * MRGcomponent ?
   * primal ou dual
   * lacanary à part ?
   * elvel of detail
   * */

  /// I have no idea what this function is supposed to do
  template<typename Int>
    void initConfigBeyerTest(LatConfig<Int>& config){}

  /// I have no idea what this function is supposed to do
  template<typename Int>
    void initConfigPalphaTest(LatConfig<Int>& config){}


} // end namespace LatticeTester

#endif

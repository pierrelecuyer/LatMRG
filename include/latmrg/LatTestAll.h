#ifndef LATTESTALL_H
#define LATTESTALL_H

#include <fnmatch.h>
#include <dirent.h>
#include <sys/types.h>
#include <cerrno>
#include <vector>
#include <string>
#include <iostream>

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/WriterRes.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/Writer.h"

#include "latmrg/LatConfig.h"
#include "latmrg/ParamReaderLat.h"
#include "latmrg/KorobovLattice.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/Formatter.h"
#include "latmrg/ReportHeaderLat.h"
#include "latmrg/ReportFooterLat.h"
#include "latmrg/ReportLat.h"
#include "latmrg/TestProjections.h"

namespace
{

  int getDir (std::string dir, std::vector <std::string> & files)
  {
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir (dir.c_str())) == NULL) {
      std::cerr << "Directory: " << dir << std::endl;
      perror ("Couldn't open the directory");
      return errno;
    }

    // Does directory name ends with /
    size_t j = dir.rfind('/');
    std::string SEP("");
    // if not, add one /
    if (dir.size() != (1 + j))
      SEP += "/";

    while ((dirp = readdir (dp)) != NULL) {
      if (0 == fnmatch("*.dat", dirp->d_name, 0))
        // keeps full name including directory name
        files.push_back (std::string (dir + SEP + dirp->d_name));
    }
    closedir (dp);
    return 0;
  }


  void eraseExtension (std::vector <std::string> & files)
  {
    for (unsigned int i = 0; i < files.size (); i++) {
      size_t j = files[i].rfind(".dat");
      if (j != std::string::npos)
        files[i].erase(j);
    }
  }


  void printFileNames (std::vector <std::string> & files)
  {
    std::cout << "----------------------------------------------" << std::endl;
    for (unsigned int i = 0; i < files.size (); i++) {
      std::cout << files[i] << std::endl;
    }
  }

}   // namespace
namespace LatMRG {

  /**
   * This class is just an auxiliary class that allows launching the spectral,
   * Beyer, or Palpha test on several different generators successively, each
   * one associated with its own data file. One can also launch these tests on all
   * the data files (with extension <tt>".dat"</tt>) in a directory. The
   * test and generator parameters in the data files must be as described in
   * program `LatMain` (see page (FIXME: page#) of this guide). In fact, the
   * `LatMain` program simply calls the methods of this class.
   */
  template<typename Int, typename Dbl>
    class LatTestAll {
      public:

        /**
         * Reads the parameters of the test and of the generator in input text file
         * `datafile`; then do the test. The data file must always have the extension
         * `".dat"`, but must be given as argument here *without extension*. For
         * example, if the data file is named `mrg.dat`, then the method must be
         * called as `doTest("mrg")`. Returns 0 if the test completed successfully;
         * returns a negative integer if there was an error.
         */
        int doTest (const char *datafile);

        /**
         * Applies the method `doTest` to all the files with extension `".dat"`
         * in directory named `dirname`. Returns 0 if all the tests completed
         * successfully; returns a non-zero integer if there was an error.
         */
        int doTestDir (const char *dirname);


        //private:
        // no longer *private* because this function *createWriter()* is sometimes
        // called outside of this class.

        /**
         * Returns a `Writer` created from the input file `infile` and the given
         * `OutputType`.
         */
        LatticeTester::Writer<Int>* createWriter (const char *infile, LatticeTester::OutputType ot);
    };

  //============================================================================

  template<typename Int, typename Dbl>
    LatticeTester::Writer<Int>* LatTestAll<Int, Dbl>::createWriter (const char *infile, LatticeTester::OutputType ot)
    {
      LatticeTester::Writer<Int> *rw = 0;
      std::string fname;

      switch (ot) {
        case LatticeTester::RES:
          fname = infile;
          fname += ".res";
          rw = new LatticeTester::WriterRes<Int> (fname.c_str ());
          break;

        case LatticeTester::TEX:
          fname = infile;
          fname += ".tex";
          // rw = new WriterTex(fname.c_str()); //EB Ne permet pas d'écrire en Tex
          break;

        case LatticeTester::TERM:
          rw = new LatticeTester::WriterRes<Int> (&std::cout);
          break;

        default:
          std::cerr << "\n*** outputType:   no such case" << std::endl;
          return 0;
      }
      return rw;
    }


  //==========================================================================

  /*
   * Reads the test parameters in infile; then do the test.
   * infile is the data file name without extension: if the data file is named
   * "poil.dat", then infile is "poil".
   * Data files must always have the extension "dat".
   */
  template<typename Int, typename Dbl>
    int LatTestAll<Int, Dbl>::doTest (const char *infile)
    {
      // Lecture des paramètres
      std::string fname (infile);
      fname += ".dat";
      ParamReaderLat<Int, Dbl> paramRdr (fname.c_str ());
      fname.clear ();

      LatConfig<Int> config;
      paramRdr.read (config);
      //config.write();

      LatticeTester::Writer<Int>* rw = createWriter (infile, config.outputType);

      LatticeTester::IntLattice<Int, Int, Dbl, Dbl> *lattice = 0;
      LatticeTester::IntLattice<Int, Int, Dbl, Dbl> *master = 0;
      LatticeTester::Lacunary<Int> *plac = 0;
      bool stationary = true;
      int toDim = config.td[1];
      int fromDim = config.td[0];
      bool memLacF = true; // Lacunary with only used lines-columns of bases
      //memLacF = false; // Lacunary with all lines-columns of bases

      if (config.J > 1) { //Several MRG
        lattice = MRGLatticeFactory<Int, Dbl>::fromCombMRG (config.comp, config.J,
            toDim, 0, config.latType, config.norm);

      } else {
        if (config.latType == PRIMEPOWER) {
          config.comp[0]->module.reduceM (config.comp[0]->a[0]);
          if (memLacF && config.lacunary)
            lattice = new MRGLatticeLac<Int, Dbl> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                config.latType, config.norm);
          else{
            lattice = new MRGLattice<Int, Dbl> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k,
                config.latType, config.norm);

          }

        } else if (config.genType[0] == MRG || config.genType[0] == LCG) {
          if (memLacF && config.lacunary){
            lattice = new MRGLatticeLac<Int, Dbl> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                config.latType, config.norm);
          }
          else{
            lattice = new MRGLattice<Int, Dbl> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k,
                config.latType, config.norm);
          }

        } else if (config.genType[0] == KOROBOV) {
          lattice = new KorobovLattice<Int, Dbl> (config.comp[0]->getM (),
              config.comp[0]->a[1], toDim, config.norm);
        } else if (config.genType[0] == RANK1) {
          stationary = false;
          lattice = new LatticeTester::Rank1Lattice<Int, Int, Dbl, Dbl> (config.comp[0]->getM (),
              config.comp[0]->a, config.comp[0]->k, config.norm);


        } else if (config.genType[0] == MMRG) {

          if (memLacF && config.lacunary) {
            lattice = new MMRGLattice<Int, Dbl> (config.comp[0]->getM(), config.comp[0]->A,
                toDim,config.comp[0]->k, config.lacunaryType,
                config.Lac, config.norm);           
          } else {
            lattice = new MMRGLattice<Int, Dbl> (config.comp[0]->getM(), config.comp[0]->A,
                toDim,config.comp[0]->k, config.norm);
          }
        }
      }

      ReportHeaderLat<Int, Dbl> header (rw, &config, lattice);
      ReportFooterLat<Int, Dbl> footer (rw);
      ReportLat<Int, Dbl> report (rw, &config, &header, &footer);

      double minVal[1 + toDim];
      LatticeTester::SetZero (minVal, toDim);

      LatticeTester::Normalizer<Dbl> *normal = 0;

      if (config.criter == LatticeTester::SPECTRAL) {
        normal = lattice->getNormalizer (config.norma, 0, config.dualF);
        // creates and returns the normalizer corresponding to config.norma
        normal->setNorm (config.norm);
      } else if (config.criter == LatticeTester::PALPHA &&
          (config.calcPalpha == LatticeTester::NORMPAL || config.calcPalpha == LatticeTester::BAL)) {
        normal = new LatticeTester::NormaPalpha<Int, Dbl> (lattice->getModulo(), config.alpha, toDim);
      }

      if (!memLacF && config.lacunary) {
        plac = new LatticeTester::Lacunary<Int> (config.Lac, toDim);
        lattice->setLac (*plac);
      }

      switch (config.criter) {
        case LatticeTester::SPECTRAL: {
                         LatTestSpectral<Int, Dbl> spectralTest (normal, lattice);
                         lattice->buildBasis (fromDim - 1);
                         spectralTest.attach (&report);

                         report.printHeader ();

                         spectralTest.setDualFlag (config.dualF);
                         spectralTest.setInvertFlag (config.invertF);
                         spectralTest.setDetailFlag (config.detailF);
                         spectralTest.setMaxAllDimFlag (true);
                         spectralTest.setMaxNodesBB (config.maxNodesBB);
                         if (1 == config.d) {
                           spectralTest.test (fromDim, toDim, minVal);
                           // lattice->write();
                           footer.setLatticeTest (&spectralTest);
                           report.printTable ();
                           report.printFooter ();

                         } else {
                           if (config.genType[0] == MRG || config.genType[0] == LCG)
                             master = new MRGLattice<Int, Dbl> (*(MRGLattice<Int, Dbl> *) lattice);
                           else if (config.genType[0] == KOROBOV)
                             master = new KorobovLattice<Int, Dbl> (*(KorobovLattice<Int, Dbl> *) lattice);
                           else if (config.genType[0] == RANK1)
                             master = new LatticeTester::Rank1Lattice<Int, Int, Dbl, Dbl> (*(LatticeTester::Rank1Lattice<Int, Int, Dbl, Dbl> *) lattice);

                           master->buildBasis (toDim);
                           TestProjections<Int, Dbl> proj (master, lattice, &spectralTest, config.td,
                               config.d);
                           proj. setOutput (rw);
                           // proj.setDualFlag (config.dualF);
                           proj.setPrintF (true);
                           double merit = proj.run (stationary, false, minVal);
                           int nbProj = proj.getNumProjections ();
                           rw->writeString ("\nMin merit:   ");
                           rw->writeDouble (sqrt (merit));
                           rw->newLine ();
                           rw->writeString ("Num projections:   ");
                           rw->writeInt (nbProj);
                           rw->newLine ();
                           // nbProj = proj.calcNumProjections(stationary, false);
                           // std::cout << "Num projections2:  " << nbProj << std::endl << std::endl;
                           delete master;
                         }
                       }
                       break;

        case LatticeTester::BEYER: {
                      LatTestBeyer<Int, Dbl> beyerTest (lattice);
                      lattice->buildBasis (fromDim - 1);
                      beyerTest.attach (&report);
                      report.printHeader ();
                      beyerTest.setDualFlag (config.dualF);
                      beyerTest.setMaxAllDimFlag (true);
                      beyerTest.setMaxNodesBB (config.maxNodesBB);
                      beyerTest.setDetailFlag (config.detailF);
                      beyerTest.test (fromDim, toDim, minVal);
                      footer.setLatticeTest (&beyerTest);
                      report.printTable ();
                      report.printFooter ();
                      //rw->writeString (lattice->toStringDualBasis ());
                    }
                    break;

        case LatticeTester::PALPHA: {
                       LatTestPalpha<Int, Dbl> palphaTest (normal, lattice);
                       palphaTest.setConfig (&config);
                       palphaTest.attach (&report);
                       report.printHeader ();
                       if (1 == config.d) {
                         palphaTest.test (fromDim, toDim, minVal);
                         footer.setLatticeTest (&palphaTest);
                         report.printTable ();
                         report.printFooter ();
                       } else {
                         MRGLattice<Int, Dbl> master = MRGLattice<Int, Dbl> (*(MRGLattice<Int, Dbl> *) lattice);
                         master.buildBasis (toDim);
                         TestProjections<Int, Dbl> proj (&master, lattice, &palphaTest, config.td,
                             config.d);
                         proj. setOutput (rw);
                         double merit = proj.run (true, false, minVal);
                         int nbProj = proj.getNumProjections ();
                         rw->writeString ("\n\nMin merit:   ");
                         rw->writeDouble (sqrt (merit));
                         rw->newLine ();
                         rw->writeString ("Num projections:   ");
                         rw->writeInt (nbProj);
                         rw->newLine ();
                         rw->newLine ();
                       }
                     }
                     break;

        default:
                     std::cerr << "Default case for config.criter" << std::endl;
                     return -1;
      }

      if (normal != 0)
        delete normal;
      if (!memLacF && config.lacunary)
        delete plac;
      delete lattice;
      delete rw;
      return 0;
    }


  //==========================================================================

  template<typename Int, typename Dbl>
    int LatTestAll<Int, Dbl>::doTestDir (const char *dirname)
    {
      std::string dir = std::string (dirname);
      std::vector <std::string> files = std::vector <std::string> ();

      getDir (dir, files);
      printFileNames (files);
      eraseExtension (files);

      int flag = 0;
      for (unsigned int i = 0; i < files.size (); i++)
        flag |= doTest (files[i].c_str());

      return flag;
    }

}
#endif

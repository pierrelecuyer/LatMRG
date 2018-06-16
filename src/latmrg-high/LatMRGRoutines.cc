// This file is part of LatMRG.
//
// LatMRG
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
#include "latmrg/LatConfig.h"
#include "latmrg/ParamReaderLat.h"
#include "latmrg/KorobovLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/LatTestAll.h"
#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/Formatter.h"
#include "latmrg/WriterRes.h"
#include "latmrg/ReportHeaderLat.h"
#include "latmrg/ReportFooterLat.h"
#include "latmrg/ReportLat.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaPalpha.h"
#include "latmrg/TestProjections.h"
#include "latmrg-high/LatMRGRoutines.h"

using namespace std;
// using namespace NTL;
using namespace LatMRG;
using namespace LatticeTester;

namespace LatMRG {

//*=============================================================================================

std::vector<double> ComputeFigureOfMerit (LatConfig<MScal>& config)
{

   //Writer* rw = createWriter (infile, config.outputType);

   std::vector<double> result;

   LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *lattice = 0;
   LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *master = 0;
   Lacunary<BScal, BVect> *plac = 0;
   int toDim = config.td[1];
   int fromDim = config.td[0];
   bool memLacF = true; // Lacunary with only used lines-columns of bases
   //memLacF = false; // Lacunary with all lines-columns of bases

   if (config.J > 1) { //Several MRG
      lattice = MRGLatticeFactory<MScal>::fromCombMRG (config.comp, config.J,
                toDim, 0, config.latType, config.norm);

   } else {
      if (config.latType == PRIMEPOWER) {
         config.comp[0]->module.reduceM (config.comp[0]->a[0]);
         if (memLacF && config.lacunary)
            lattice = new MRGLatticeLac<MScal> (config.comp[0]->module.mRed,
               config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                                         config.latType, config.norm);
         else{
            lattice = new MRGLattice<MScal> (config.comp[0]->module.mRed,
               config.comp[0]->a, toDim, config.comp[0]->k,
                                      config.latType, config.norm);

         }

      } else if (config.genType[0] == MRG || config.genType[0] == LCG) {
     
         if (memLacF && config.lacunary){            
            lattice = new MRGLatticeLac<MScal> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k, config.Lac,
                config.latType, config.norm);
         }
         else{
            lattice = new MRGLattice<MScal> (config.comp[0]->module.mRed,
                config.comp[0]->a, toDim, config.comp[0]->k,
                config.latType, config.norm);
         }

      } else if (config.genType[0] == KOROBOV) {
         lattice = new KorobovLattice<MScal> (config.comp[0]->getM (),
            config.comp[0]->a[1], toDim, config.norm);
      } else if (config.genType[0] == RANK1) {
         lattice = new Rank1Lattice<MScal, MVect, BScal, BVect, BMat, NScal, NVect, RScal> (config.comp[0]->getM (),
            config.comp[0]->a, config.comp[0]->k, config.norm);
      

      } else if (config.genType[0] == MMRG) {
         if (memLacF && config.lacunary) {
            lattice = new MMRGLattice<MScal> (config.comp[0]->getM(), config.comp[0]->A,
                             toDim,config.comp[0]->k, config.lacunaryType,
                             config.Lac, config.norm);           
         } else {
            lattice = new MMRGLattice<MScal> (config.comp[0]->getM(), config.comp[0]->A,
                             toDim,config.comp[0]->k, config.norm);
         }
      }
   }

   double minVal[1 + toDim];
   SetZero (minVal, toDim);

   Normalizer<RScal> *normal = 0;

   if (config.criter == SPECTRAL) {
      normal = lattice->getNormalizer (config.norma, 0, config.dualF);
      // creates and returns the normalizer corresponding to config.norma
      normal->setNorm (config.norm);
   } else if (config.criter == PALPHA &&
              (config.calcPalpha == NORMPAL || config.calcPalpha == BAL)) {
      normal = new NormaPalpha<MScal, RScal> (lattice->getModulo(), config.alpha, toDim);
   }

   if (!memLacF && config.lacunary) {
      plac = new Lacunary<BScal, BVect> (config.Lac, toDim);
      lattice->setLac (*plac);
   }

   switch (config.criter) {
   case SPECTRAL: {
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
            // lattice->write();
            //footer.setLatticeTest (&spectralTest);
            //report.printTable ();
            //report.printFooter ();
         } else {
            if (config.genType[0] == MRG || config.genType[0] == LCG)
               master = new MRGLattice<MScal> (*(MRGLattice<MScal> *) lattice);
            else if (config.genType[0] == KOROBOV)
               master = new KorobovLattice<MScal> (*(KorobovLattice<MScal> *) lattice);
            else if (config.genType[0] == RANK1)
               master = new Rank1Lattice<MScal, MVect, BScal, BVect, BMat, NScal, NVect, RScal> (*(Rank1Lattice<MScal, MVect, BScal, BVect, BMat, NScal, NVect, RScal> *) lattice);

            master->buildBasis (toDim);
            TestProjections proj (master, lattice, &spectralTest, config.td,
                                  config.d);
            //proj. setOutput (rw);
            proj.setDualFlag (config.dualF);
            proj.setPrintF (true);
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

   case BEYER: {
         LatTestBeyer<MScal, BScal, BVect, BMat, NScal, NVect, RScal> beyerTest (lattice);
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

   case PALPHA: {
         LatTestPalpha<MScal> palphaTest (normal, lattice);
         palphaTest.setConfig (&config);
         //palphaTest.attach (&report);
         //report.printHeader ();
         if (1 == config.d) {
            palphaTest.test (fromDim, toDim, minVal);
            //footer.setLatticeTest (&palphaTest);
            //report.printTable ();
            //report.printFooter ();
         } else {
            MRGLattice<MScal> master = MRGLattice<MScal> (*(MRGLattice<MScal> *) lattice);
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
      cerr << "Default case for config.criter" << endl;
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

//*=============================================================================================

void printResult(const std::vector<double> & result, const int & fromDim)
{
   cout << "Result: " << endl;
   for (unsigned int i = 0; i < result.size(); ++i)
      cout << "  dim " << fromDim+i << " = " << sqrt(result[i]) << endl;

}

//*=============================================================================================

void initConfigSpectralTest(LatConfig<MScal>& config)
{}

//*=============================================================================================

void initConfigBeyerTest(LatConfig<MScal>& config)
{}

//*=============================================================================================

void initConfigPalphaTest(LatConfig<MScal>& config)
{}

//*=============================================================================================

}

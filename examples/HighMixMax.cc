/*
 * This file showcases the usage of the LatMRG API at a high level. It contains
 * the computations of the MIXMAX article (to be published). This program does
 * something very similar to LowMixMax.cc but has hard coded parameters instead
 * of taking input from standard input.
 *
 * To compile this program you should use something similar to :
 *      g++ -I./include -I./latticetester/include -c \
 *      ./examples/HighMixMax.cc -L./lib -llatmrg -llatticetester \
 *      -L/usr/local/lib -lntl -lgmp -o ./bin/examples/HighMixMax
 * with LatMRG as working directory with the LatMRG and LatticeTester library 
 * built.
 *
 * You can get the tests in the MixMax article by calling :
 *      ./bin/examples/HighMixMax
 * in the LatMRG directory (assuming you compiled with the command above). 
 * Since the parameters are hard coded, the output is self explanatory if you
 * peak at the article.
 *
 * Authors: - Paul Wambergue
 *          - Marc-Antoine Savard
 * */

// We need to handle large integers, but we don't need to be that precise when 
// we do floating point
#define NTL_TYPES_CODE 2

// Include Header
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

// Include LatMRG Header
#include "latmrg/LatConfig.h"
#include "latmrg-high/LatMRGRoutines.h"
#include "latmrg/mrgtypes/MixmaxMMRG.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;
using namespace LatMRG;


void doTests(MScal m, int k, MScal d, MScal c) {
   MixmaxMMRG mixmax (m, k, d, c);
   //cout << "A = \n" << mixmax.getMatrix() << endl;

   LatConfig<MScal> latconfig;

   latconfig.J = 1;
   latconfig.setJ(latconfig.J);
   for (int i = 0; i < latconfig.J; i++) {
      latconfig.genType[i] = MMRG; //PW_TODO moche
      latconfig.comp[i] = new MRGComponent<MScal>(m, mixmax.getMatrix(), k);
   }
   
   //latconfig.comp[0] = new MRGComponent (m, A, k);
   //latconfig.genType[0] = MMRG;

   latconfig.d = 1;
   latconfig.td = new int[latconfig.d];

   latconfig.criter = SPECTRAL;
   latconfig.norma = BESTLAT;
   latconfig.norm = L2NORM;
   latconfig.dualF = true;
   latconfig.invertF = false;
   latconfig.latType = FULL;
   latconfig.lacunary = true;
   latconfig.lacunaryType = ARBITRARYINDICES;
   latconfig.numberLacIndices = 7;
   latconfig.maxNodesBB = 1000000;

   // Proposition 4
   latconfig.td[0] = 5;
   latconfig.td[1] = 5;
   latconfig.Lac.resize(5);
   latconfig.Lac[0]=4;
   latconfig.Lac[1]=5;
   latconfig.Lac[2]=k+3;
   latconfig.Lac[3]=k+4;
   latconfig.Lac[4]=k+5;

   std::vector<double> FoM;
   FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);

   printResult(FoM, latconfig.td[0]);

   // Proposition 6
   latconfig.td[0] = 7;
   latconfig.td[1] = 7;
   latconfig.Lac.resize(7);
   latconfig.Lac[0]=4;
   latconfig.Lac[1]=5;
   latconfig.Lac[2]=6;
   latconfig.Lac[3]=k+3;
   latconfig.Lac[4]=k+4;
   latconfig.Lac[5]=k+5;
   latconfig.Lac[6]=k+6;
   
   FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);

   printResult(FoM, latconfig.td[0]);
}
//==================================================================================


int main ()
{
   // four-family parameters DIMENSION 8
   MScal m = power_ZZ(2,61) - 1; // p
   int k = 8; // N
   MScal d (0); // s
   MScal c = power_ZZ(2,53) + 1; //m
   doTests(m, k, d, c);
   
   // four-family parameters DIMENSION 240
   m = power_ZZ(2,61) - 1; // p
   k = 240; // N
   d = NTL::ZZ(487013230256099140); // s
   c = power_ZZ(2,51) + 1; //m
   doTests(m, k, d, c);
   //
   // four-family parameters DIMENSION 240
   m = power_ZZ(2,61) - 1; // p
   k = 17; // N
   d = NTL::ZZ(0); // s
   c = power_ZZ(2,36) + 1; //m
   doTests(m, k, d, c);
   

   return 0;
}



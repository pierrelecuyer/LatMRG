//
//  routinesTest.cpp
//  main program to test main functions for
//  high level usage of LatticeTester
//

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
#include "latmrg/LatMRGRoutines.h"
#include "latmrg/LatConfig.h"
#include "latmrg/MixmaxMMRG.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;
using namespace LatMRG;


//==================================================================================


int main ()
{
   // four-family parameters DIMENSION 8
   /*
   MScal m = power_ZZ(2,61) - 1; // p
   int k = 8; // N
   MScal d (0); // s
   MScal c = power_ZZ(2,53) + 1; //m
   */
   
   // four-family parameters DIMENSION 240
   MScal m = power_ZZ(2,61) - 1; // p
   int k = 240; // N
   MScal d (487013230256099140); // s
   MScal c = power_ZZ(2,51) + 1; //m
   
   MixmaxMMRG mixmax (m, k, d, c);
   cout << "A = \n" << mixmax.getMatrix() << endl;

   LatConfig latconfig;

   latconfig.J = 1;
   latconfig.setJ(latconfig.J);
   for (int i = 0; i < latconfig.J; i++) {
      latconfig.genType[i] = MMRG; //PW_TODO moche
      latconfig.comp[i] = new MRGComponent(m, mixmax.getMatrix(), k);
   }
   
   //latconfig.comp[0] = new MRGComponent (m, A, k);
   //latconfig.genType[0] = MMRG;

   latconfig.d = 1;
   latconfig.td = new int[latconfig.d];
   latconfig.td[0] = 7;
   latconfig.td[1] = 7;

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

   latconfig.Lac.resize(latconfig.numberLacIndices);
   // DIMENSION 8 - first part
   /*
   latconfig.Lac[0]=4;
   latconfig.Lac[1]=5;
   latconfig.Lac[2]=6;
   latconfig.Lac[3]=11;
   latconfig.Lac[4]=12;
   latconfig.Lac[5]=13;
   latconfig.Lac[6]=14;
   */
   
   // DIMENSION 8 - second part
   /*
   latconfig.Lac[0]=4;
   latconfig.Lac[1]=5;
   latconfig.Lac[2]=11;
   latconfig.Lac[3]=12;
   latconfig.Lac[4]=13;
   */
   
   // DIMENSION 240
   latconfig.Lac[0]=4;
   latconfig.Lac[1]=5;
   latconfig.Lac[2]=6;
   latconfig.Lac[3]=243;
   latconfig.Lac[4]=244;
   latconfig.Lac[5]=245;
   latconfig.Lac[6]=246;
   
   std::vector<double> FoM;
   FoM = ComputeFigureOfMerit(latconfig);

   printResult(FoM, latconfig.td[0]);

   return 0;
}



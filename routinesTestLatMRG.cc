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

using namespace std;
using namespace NTL;
using namespace LatticeTester;
using namespace LatMRG;


//==================================================================================


int main ()
{

   MScal m = power_ZZ(2,61) - 1;
   int k = 8;
   MMat A;
   A.resize(8, 8);

   A[0][0]=1; A[0][1]=1;                 A[0][2]=1;                 A[0][3]=1;                 A[0][4]=1;                 A[0][5]=1;                 A[0][6]=1;                A[0][7]=1;
   A[1][0]=1; A[1][1]=2;                 A[1][2]=1;                 A[1][3]=1;                 A[1][4]=1;                 A[1][5]=1;                 A[1][6]=1;                A[1][7]=1;
   A[2][0]=1; A[2][1]=9007199254740995;  A[2][2]=2;                 A[2][3]=1;                 A[2][4]=1;                 A[2][5]=1;                 A[2][6]=1;                A[2][7]=1;
   A[3][0]=1; A[3][1]=18014398509481988; A[3][2]=9007199254740995;  A[3][3]=2;                 A[3][4]=1;                 A[3][5]=1;                 A[3][6]=1;                A[3][7]=1;
   A[4][0]=1; A[4][1]=27021597764222981; A[4][2]=18014398509481988; A[4][3]=9007199254740995;  A[4][4]=2;                 A[4][5]=1;                 A[4][6]=1;                A[4][7]=1;
   A[5][0]=1; A[5][1]=36028797018963974; A[5][2]=27021597764222981; A[5][3]=18014398509481988; A[5][4]=9007199254740995;  A[5][5]=2;                 A[5][6]=1;                A[5][7]=1;
   A[6][0]=1; A[6][1]=45035996273704967; A[6][2]=36028797018963974; A[6][3]=27021597764222981; A[6][4]=18014398509481988; A[6][5]=9007199254740995;  A[6][6]=2;                A[6][7]=1;
   A[7][0]=1; A[7][1]=54043195528445960; A[7][2]=45035996273704967; A[7][3]=36028797018963974; A[7][4]=27021597764222981; A[7][5]=18014398509481988; A[7][6]=9007199254740995; A[7][7]=2;

   cout << "A = \n" << A << endl;

   LatConfig latconfig;


   latconfig.J = 1;
   latconfig.setJ(latconfig.J);
   for (int i = 0; i < latconfig.J; i++) {
      latconfig.genType[i] = MMRG; //PW_TODO moche
      latconfig.comp[i] = new MRGComponent(m, A, k);
   }
   
   //latconfig.comp[0] = new MRGComponent (m, A, k);
   //latconfig.genType[0] = MMRG;

   latconfig.d = 1;
   latconfig.td = new int[latconfig.d];
   latconfig.td[0] = 7;
   latconfig.td[1] = 10;

   latconfig.criter = SPECTRAL;
   latconfig.norma = BESTLAT;
   latconfig.norm = L2NORM;
   latconfig.dualF = true;
   latconfig.invertF = false;
   latconfig.latType = FULL;
   latconfig.lacunary = true;
   latconfig.lacunaryType = ARBITRARYINDICES;
   latconfig.numberLacIndices = 10;
   latconfig.maxNodesBB = 1000000;

   latconfig.Lac.resize(latconfig.numberLacIndices);
   latconfig.Lac[0]=2;
   latconfig.Lac[1]=3;
   latconfig.Lac[2]=4;
   latconfig.Lac[3]=5;
   latconfig.Lac[4]=6;
   latconfig.Lac[5]=7;
   latconfig.Lac[6]=10;
   latconfig.Lac[7]=11;
   latconfig.Lac[8]=12;
   latconfig.Lac[9]=13;

   std::vector<double> FoM;
   FoM = applyTest(latconfig);

   printResult(FoM, latconfig.td[0]);

   return 0;
}



//
//  main_for_mmrg.cc
//  LatMRG
//
//  Main program to test the new MMRG class
//


// select pre compiling options
//----------------------------------------------------------------------------------------

#define PRINT_CONSOLE

//----------------------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <num.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

#include "latticetester/Types.h"
#include "latticetester/Util.h"

#include "latmrg/IntLattice.h"
#include "latmrg/LatTestAll.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MMRGLattice.h"

//#include "SimpleMRG.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;
using namespace LatMRG;

//*=======================================================================================

int main ()
{
   
   // modulus
   //ZZ m = power_ZZ(2,5)-1; //=31
   ZZ m = power_ZZ(2,7)-1; //=127
   
   // dimension of projection
   int dimension = 7;
   
   // generator matrix
   int r = 4;
   MMat A;
   A.SetDims(r, r);
   A[0][0]=2; A[0][1]=0; A[0][2]=2; A[0][3]=9;
   A[1][0]=3; A[1][1]=1; A[1][2]=0; A[1][3]=2;
   A[2][0]=4; A[2][1]=5; A[2][2]=1; A[2][3]=3;
   A[3][0]=7; A[3][1]=3; A[3][2]=5; A[3][3]=0;
   
   cout << "det(A) = " << determinant(A) << endl;
   
   // MMRG
   MMRGLattice myMMRG (m, A, 3*r, r, FULL, L2NORM);
   cout << "A = \n" << myMMRG.toStringGeneratorMatrix() << endl;
   
   myMMRG.buildBasis(dimension);
   cout << "Primal basis = \n" << myMMRG.getBasis() << endl;
   cout << "Dual basis = \n" << myMMRG.getDualBasis() << endl;
   
   cout << "V*transpose(W) = \n" << myMMRG.getBasis() * transpose(myMMRG.getDualBasis()) << endl;
   
   
   myMMRG.incrementDimBasis();
   cout << "\n----- Increment basis dimension -----" << endl;
   cout << "Primal basis = \n" << myMMRG.getBasis() << endl;
   cout << "Dual basis = \n" << myMMRG.getDualBasis() << endl;
   
   
   return 0;
}


//*=======================================================================================

/*
 
 int main (int argc, char *argv[])
 {
 // Erwan
 //string testLocation = "/Users/Erwan1/projects/github/LatMRG";
 //string testLocation = "/Users/Erwan1/projects/github/LatMRG/latZZDD_test1";
 
 // Paul
 //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles";
 string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/latZZDD_test1";
 
 struct stat buf; //properties of a file or directory
 LatMRG::LatTestAll testall;
 int status = 0;
 
 stat(testLocation.c_str(), &buf);
 
 if (0 != S_ISDIR(buf.st_mode)) //directory
 status |= testall.doTestDir (testLocation.c_str());
 else { //file
 string dataname(testLocation.c_str());
 dataname.append(".dat");
 stat(dataname.c_str(), &buf);
 
 status |= testall.doTest (testLocation.c_str());
 }
 return 0;
 }
 
 */
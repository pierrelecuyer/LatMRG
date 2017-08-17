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

void SavvidyMatrix(mat_ZZ& A, int N, int s, int m, int b)
{
   A.kill();
   A.SetDims(N,N);
   for (int j = 1; j < N; j ++) {
      for (int i = j+1; i < N; i++)
         A[i][j] = (i-j+2) * m + b;
   }
   for (int i = 0; i < N; i ++)
      A[i][0] = 1;
   for (int i = 1; i < N; i++)
      A[i][i] = 2;
   for (int i = 0; i < N; i ++) {
      for (int j = i+1; j < N; j ++)
         A[i][j] = 1;
   }
   A[2][1] += s;
}

void LacunaryMatrix (mat_ZZ& B, mat_ZZ& A, vec_ZZ& lacunaryIndices)
{
   int L = (int) lacunaryIndices.length();
   int N = (int) A.NumCols();

   B.SetDims(L, N);
   for (int k = 0; k < L; k++)
      B[k][conv<int>(lacunaryIndices[k])-1] = 1;
}

void printMatrixForInputFile (mat_ZZ A)
{

    for (int i = 0; i < A.NumRows(); i++) {
        for (int j = 0; j < A.NumCols(); j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

//*=======================================================================================

#if 0
int main ()
{
   
   // modulus
   //ZZ m = power_ZZ(2,3)-1; 
   ZZ m = power_ZZ(2,7)-1; //=127
   
   // dimension of projection
   int dimension = 7;
   
   // generator matrix
   int r = 5;
   MMat A;
   SavvidyMatrix(A, r, -1, 1, 0);

   // MMRG
   MMRGLattice myMMRG (m, A, dimension, r, L2NORM, FULL);
   cout << "A = \n" << myMMRG.toStringGeneratorMatrix() << endl;

   // MMRG lacunary
   vec_ZZ lacunaryIndices;
   lacunaryIndices.SetLength(3);
   lacunaryIndices[0] = 1;
   lacunaryIndices[1] = 3;
   lacunaryIndices[2] = 5;
   //lacunaryIndices[3] = 4;
   //lacunaryIndices[4] = 5;
   MMat B;
   LacunaryMatrix(B, A, lacunaryIndices);
   cout << "lac B = \n" << B << endl;
   MMRGLattice myLacMMRG (m, A, dimension, r, L2NORM, FULL);
   cout << "lac A = \n" << myLacMMRG.toStringGeneratorMatrix() << endl;



   cout << "\n----- Build basis non lac -----" << endl;
   myMMRG.buildNonLacunaryBasis(dimension);
   cout << "Primal basis = \n"; 
   printMatrixForInputFile(myMMRG.getBasis());
   cout << "Dual basis = \n"; 
   printMatrixForInputFile(myMMRG.getDualBasis());
   cout << "m_dualVecNorm = " << myMMRG.getDualVecNorm() << endl;
   //cout << "V*transpose(W) = \n" << myMMRG.getBasis() * transpose(myMMRG.getDualBasis()) << endl;

   cout << "\n----- Build basis lac -----" << endl;
   myLacMMRG.buildLacunaryBasis(dimension, B);
   cout << "Primal basis = \n";
   printMatrixForInputFile(myLacMMRG.getBasis());
   cout << "Dual basis = \n"; 
   printMatrixForInputFile(myLacMMRG.getDualBasis());
   cout << "m_dualVecNorm = " << myLacMMRG.getDualVecNorm() << endl;
   //cout << "V*transpose(W) = \n" << myLacMMRG.getBasis() * transpose(myLacMMRG.getDualBasis()) << endl;



   cout << "\n----- Increment dim non lac -----" << endl;
   for (int j = 0; j < 10; j++)
      myMMRG.incrementDimBasis();
   cout << "Primal basis = \n"; 
   printMatrixForInputFile(myMMRG.getBasis());
   //cout << "Dual basis = \n" << myMMRG.getDualBasis() << endl;
   //cout << "V*transpose(W) = \n" << myMMRG.getBasis() * transpose(myMMRG.getDualBasis()) << endl;
   
   cout << "\n--- Increment dim lac ---" << endl;
   for (int j = 0; j < 10; j++)
      myLacMMRG.incrementDimLacunaryBasis(B);
   cout << "Primal basis = \n";
   printMatrixForInputFile(myLacMMRG.getBasis());
   //cout << "Dual basis = \n" << myLacMMRG.getDualBasis() << endl;
   //cout << "V*transpose(W) = \n" << myLacMMRG.getBasis() * transpose(myLacMMRG.getDualBasis()) << endl;
   
   return 0;
}
#endif

//*=======================================================================================

#if 0
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
 
#endif
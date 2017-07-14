//
// Brouillon de travail
// main program to build lattice basis for Matrix-MRG
// 

#include <iostream>
#include <map>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"
#include "latticetester/Types.h"

#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

using namespace std;
using namespace NTL;
using namespace LatticeTester;


void SavvidyMatrix(mat_ZZ& A, int N, int s, int m, int b) {

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


mat_ZZ Dualize (const mat_ZZ V, const ZZ modulus, const int k)
{
   mat_ZZ W;
   W.SetDims(V.NumRows(), V.NumRows());

   transpose(W,-V);
   long rmax = k;
   if(k>V.NumRows()){ rmax = V.NumRows(); }
   for (int i = 0; i < rmax; i++)
      W[i][i] = modulus;
   for (int i = k; i < V.NumRows(); i++)
      W[i][i] = 1;

   return W;
}


mat_ZZ buildBasis (mat_ZZ& A, int dimension, ZZ modulus)
{
    BMat B;
    int sizeA = A.NumCols();
    B.resize(dimension, dimension);

    // filling in the diagonal of B
    for (int i = 0; i < sizeA; i++)
        B[i][i] = 1;
    for (int i = sizeA; i < dimension; i++)
        B[i][i] = modulus;

    // filling in the values specific to the matrix A
    ZZ_p::init(modulus);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

    int maxIter = dimension/sizeA;

    for (int k = 1; k < maxIter+1; k++) {
        // calcul de transpose(A^k)
        temp *= conv<mat_ZZ_p>(transpose(A)); 

        if (k == maxIter) { //on complète le bout de la matrice B
            int residu = dimension - maxIter * sizeA;
            for (int i = 0; i < sizeA; i++) {
                for (int j = 0; j < residu; j ++)
                    B[i][k*sizeA +j] = conv<ZZ>(temp[i][j]);
            }
        } else {
            for (int i = 0; i < sizeA; i++) {
                for (int j = 0; j < sizeA; j ++)
                    B[i][k*sizeA +j] = conv<ZZ>(temp[i][j]);
            }
        }
    }

    return B;
}

void incrementDimension (mat_ZZ& B, mat_ZZ& A, ZZ modulus)
{

    int oldDimension = B.NumRows();
    int newDimension = oldDimension+1;
    int sizeA = A.NumRows();

    mat_ZZ Btemp = B;
    B.SetDims(oldDimension+1, oldDimension+1);

    for (int i = 0; i < oldDimension; i++) {
        for (int j = 0; j < oldDimension; j++)
            B[i][j] = Btemp[i][j];
    }

    // calcul de la puissance a A en cours
    int n = floor(newDimension / sizeA);
    ZZ_p::init(modulus);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++) {
        for (int j = 0; j < sizeA; j++)
            temp[i][j] = conv<ZZ_p>(B[i][j]);
    }

    // étape couteuse qui pourrait etre raccourcie en stockant A^k
    //--------------------------------------------------------------------
    for (int k = 1; k < n+1; k++) {
        temp *= conv<mat_ZZ_p>(transpose(A)); 
    }
    //--------------------------------------------------------------------

    // mise à jour de la nouvelle colonne de B
    for (int i = 0; i < sizeA; i++) {
        B[i][newDimension-1] = conv<ZZ>(temp[i][newDimension - n*sizeA -1]);
    }

    B[newDimension-1][newDimension-1] = modulus;

}


mat_ZZ buildBasisMultiple (mat_ZZ& A, int multiple, ZZ modulus)
{
    BMat B;
    int sizeA = A.NumCols();
    B.resize(multiple*sizeA, multiple*sizeA);

    // filling in the diagonal of B
    for (int i = 0; i < sizeA; i++)
        B[i][i] = 1;
    for (int i = sizeA; i < multiple*sizeA; i++)
        B[i][i] = modulus;

    // filling in the values specific to the matrix A
    ZZ_p::init(modulus);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

    for (int k = 1; k < multiple; k++) {
        // calcul de transpose(A^k)
        temp *= conv<mat_ZZ_p>(transpose(A)); 

        for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < sizeA; j ++)
                B[i][k*sizeA +j] = conv<ZZ>(temp[i][j]);
        }
    }

    return B;
}
    

void printMatrixForInputFile (mat_ZZ & A)
{

    for (int i = 0; i < A.NumRows(); i++) {
        for (int j = 0; j < A.NumCols(); j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}


//*===========================================================*

int main () 
{

    cout << "\nTest Savvidy\n" << endl;
    mat_ZZ Atest;
    SavvidyMatrix(Atest, 20, 0, 1, 0);
    cout << "Atest =\n" << Atest << endl;
    cout << "\n----------------------\n" << endl;



    // dimension parameters
    int r = 10;
    int n = 4;

    // modulus
    //ZZ p = conv<ZZ>(101);
    ZZ p = power_ZZ(2,7) - 1; //=127

    /*
    // matrix fir usual MRG

    // coefficients (a_i)pyht
    BVect a;
    a.resize(r);
    a[0]=2;
    a[1]=5;
    a[2]=17;
    a[3]=3;
    a[4]=34;

    BMat A;
    A.resize(r,r);
    for (int i = 0; i < r-1; i++)
        A[i][i+1] = 1;
    for (int j = 0; j < r; j++)
        A[r-1][j] = a[r-1-j];
    cout << "A = \n" << A << endl;
    */

    // savvidy matrix 
    mat_ZZ A;
    SavvidyMatrix(A, r, -1, 1, 0);
    cout << "A = \n" << A << endl;

    cout << "tr(A) = \n" << transpose(A) << endl;
    cout << "tr(A^2) = \n" << transpose(A*A) << endl;

    // matrix containing a basis of the lattice
    BMat B;
    B.resize(n*r,n*r);
    for (int i = 0; i < r; i++)
        B[i][i] = 1;
    for (int i = r; i < n*r; i++)
        B[i][i] = p;
    ZZ_p::init(p);
    mat_ZZ_p temp;
    temp.SetDims(r,r);
    for (int i = 0; i < r; i++)
        temp[i][i] = 1;

    for (int k = 1; k < n; k++) {
        // calcul de transpose(A^k)
        temp *= conv<mat_ZZ_p>(transpose(A)); 

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < r; j ++)
                B[i][k*r+j] = conv<ZZ>(temp[i][j]);
        }
    }

    printMatrixForInputFile(B);


    //cout << "\nB = \n" << B << endl;
    //mat_ZZ Bbis = buildBasisMultiple(A, n, p);
    //cout << "\nBbis = \n" << Bbis << endl;

    /*
    for (int i = r; i < n*r; i++) {
        mat_ZZ Btemp;
        Btemp = buildBasis(A, i, p);
        cout << "\nB(" << i << ") = \n" << Btemp << endl;
    }
    */


    cout << "Avec r = 9 pas de soucis, mais avec r =10 (la dimension 30 en est donc un multiple) ça pose problème. ";
    cout << "Ça vient de la partie *** calcul de la puissance de A en cours***" << endl;

    mat_ZZ Bbis;
    Bbis = buildBasis(A, 30, p);
    mat_ZZ B_inc;
    B_inc = buildBasis(A, 29, p);
    //cout << "\nB_inc before = \n" << B_inc << endl;
    incrementDimension(B_inc, A, p);
    //cout << "\nB = \n" << Bbis << endl;
    //cout << "\nB_inc after = \n" << B_inc << endl;

    /*
    mat_ZZ Bdual;
    Bdual = Dualize(B, p, r);
    cout << "Test dual = \n" << B*transpose(Bdual) << endl;

    cout << "Bdual = \n" << Bdual << endl;
    printMatrixForInputFile(Bdual);
    */

    return 0;
}





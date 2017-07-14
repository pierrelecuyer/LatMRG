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


mat_ZZ buildBasis (mat_ZZ& A, int dimension, ZZ modulus, int factor)
{
    BMat B;
    int sizeA = A.NumCols();
    B.resize(dimension, dimension);

    // filling in the diagonal of B
    for (int i = 0; i < sizeA; i++)
        B[i][i] = factor;
    for (int i = sizeA; i < dimension; i++)
        B[i][i] = modulus;

    // filling in the values specific to the matrix A
    ZZ_p::init(modulus);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++)
        temp[i][i] = factor;

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
    B.SetDims(newDimension, newDimension);

    for (int i = 0; i < oldDimension; i++) {
        for (int j = 0; j < oldDimension; j++)
            B[i][j] = Btemp[i][j];
    }

    // calcul de la puissance a A en cours
    int n = floor(oldDimension / sizeA);
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

vec_ZZ getSubLine(mat_ZZ& B, int lign, int jMin, int jMax)
{
    // both jMin and jMax are included
    vec_ZZ vec;
    vec.SetLength(jMax-jMin+1);
    for (int i = 0; i < (jMax-jMin+1); i++)
        vec[i] = B[lign][jMin+i];
    return vec;
}

void incrementDimensionTest (mat_ZZ& B, mat_ZZ& A, ZZ modulus)
{
    int oldDimension = B.NumRows();
    int newDimension = oldDimension+1;
    int sizeA = A.NumRows();

    mat_ZZ Btemp = B;
    B.SetDims(newDimension, newDimension);

    for (int i = 0; i < oldDimension; i++) {
        for (int j = 0; j < oldDimension; j++)
            B[i][j] = Btemp[i][j];
    }

    // calcul de la puissance a A en cours
    int n = floor(oldDimension / sizeA);
    ZZ_p::init(modulus);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++) {
            temp[i][i] = 1;
    }

    // étape couteuse qui pourrait etre raccourcie en stockant A^k
    //--------------------------------------------------------------------
    for (int k = 1; k < n+1; k++) {
        temp *= conv<mat_ZZ_p>(transpose(A)); 
    }
    //--------------------------------------------------------------------

    // mise à jour de la nouvelle colonne de B
    vec_ZZ initialState;
    for (int i = 0; i < oldDimension; i++) {
        initialState = getSubLine(B, i, 0, sizeA-1);
        initialState = conv<vec_ZZ>( transpose(temp) * conv<vec_ZZ_p>(initialState) );
        B[i][newDimension-1] = initialState[newDimension - n*sizeA -1];
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
    

void printMatrixForInputFile (mat_ZZ A)
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

    // dimension parameters
    int r = 3;
    int dimension = 12;

    // modulus
    //ZZ p = conv<ZZ>(101);
    ZZ p = power_ZZ(2,7) - 1; //=127

    /*
    // matrix for usual MRG
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




    /*
    // savvidy matrix 
    mat_ZZ A;
    SavvidyMatrix(A, r, -1, 1, 0);
    //cout << "A = \n" << A << endl;

    mat_ZZ B;
    B = buildBasis(A, dimension+1, p, 1);

    mat_ZZ B1;
    B1 = buildBasis(A, dimension, p, 7);
    cout << "B1 = \n" << B1 << endl;
    incrementDimension (B1, A, p);
    cout << "B1 = +1 dim\n" << B1 << endl;

    mat_ZZ B2;
    B2 = buildBasis(A, dimension, p, 6);
    cout << "B2 = \n" << B2 << endl;
    incrementDimension (B2, A, p);
    cout << "B2 = +1 dim\n" << B2 << endl;

    printMatrixForInputFile(B-(B1-B2));


    cout << "------------------------------------------------" << endl;

    mat_ZZ Bsure;
    Bsure = buildBasis(A, dimension, p, 2);
    cout << "Bsure = \n" << Bsure << endl;
    mat_ZZ Btest;
    Btest = buildBasis(A, dimension, p, 2);
    cout << "Btest = \n" << Btest << endl;

    cout << "B = \n" << B << endl;

    incrementDimension(Bsure, A, p);
    cout << "Bsure = \n" << Bsure << endl;

    incrementDimensionTest(Btest, A, p);
    cout << "Btest = \n" << Btest << endl;

    cout << "------------------------------------------------" << endl;
    */




    
    // savvidy matrix 
    mat_ZZ A;
    SavvidyMatrix(A, r, -1, 1, 0);
    cout << "A = \n" << A << endl;

    mat_ZZ init;
    init.SetDims(r,r);
    init[0][0]=2; init[0][1]=2; init[0][2]=0;
    init[1][0]=8; init[1][1]=9; init[1][2]=12;
    init[2][0]=3; init[2][1]=0; init[2][2]=2;
    cout << "det(init) = " << determinant(init) << endl;
    cout << "init = \n" << init << endl;

    cout << "init * tr(A) = \n" << init * transpose(A) << endl;
    cout << "init * tr(A^2) = \n" << init * transpose(A*A) << endl;
    cout << "init * tr(A^3) = \n" << init * transpose(A*A*A) << endl;

    mat_ZZ B;
    B.SetDims(r+1,r+1);
    B[0][0]=2; B[0][1]=2; B[0][2]=0; B[0][3]=15;
    B[1][0]=8; B[1][1]=9; B[1][2]=12; B[1][3]=1;
    B[2][0]=3; B[2][1]=0; B[2][2]=2; B[2][3]=16;
    B[3][0]=2; B[3][1]=4; B[3][2]=8; B[3][3]=10;
    cout << "det(B) = " << determinant(B) << endl;
    cout << "B = \n" << B << endl;

    mat_ZZ Bbis = B;
    cout << "Bbis = \n" << Bbis << endl;

    cout << "\n Incrementing dimension:" << endl;

    for (int k = 0; k < 10; k++) {

        incrementDimension(Bbis, A, p);
        cout << "Bbis = \n" << Bbis << endl;

        incrementDimensionTest(B, A, p);
        cout << "B = \n" << B << endl;

        cout << "---------------------" << endl;

    }
    


    /*
    mat_ZZ B;
    B = buildBasis(A, dimension, p);
    mat_ZZ B_inc;
    B_inc = buildBasis(A, dimension-1, p);
    cout << "\nB_inc before = \n" << B_inc << endl;
    incrementDimension(B_inc, A, p);
    cout << "\nB = \n" << B << endl;
    cout << "\nB_inc after = \n" << B_inc << endl;
    */

    //printMatrixForInputFile(B_inc - B);


    /*
    mat_ZZ Bdual;
    Bdual = Dualize(B, p, r);
    cout << "Test dual = \n" << B*transpose(Bdual) << endl;

    cout << "Bdual = \n" << Bdual << endl;
    printMatrixForInputFile(Bdual);
    */

    return 0;
}





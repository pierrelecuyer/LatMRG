#include <iostream>
#include "latmrg/LatTestSpectral.h"

using namespace std;
using namespace LatMRG;

void SavvidyMatrix(mat_ZZ_p& A, ZZ p, int N, ZZ s, ZZ m, ZZ b) {
	ZZ_p::init(p);
    A.kill();
    A.SetDims(N,N);
    for (int j = 1; j < N; j ++) {
        for (int i = j+1; i < N; i++)
            A[i][j] = conv<ZZ_p>((i-j+2) * m + b);
    }
    for (int i = 0; i < N; i ++)
        A[i][0] = 1;
    for (int i = 1; i < N; i++)
        A[i][i] = 2;
    for (int i = 0; i < N; i ++) {
        for (int j = i+1; j < N; j ++)
            A[i][j] = 1;
    }
    A[2][1] += conv<ZZ_p>(s);
}

//==========================================================================

int main ()
{
	// parameters of the test
	CriterionType criterion = BEYER;		
	NormType norm = L2NORM;
	bool ReadGenFile = false;						
	int numberOfComponents = 1;
	GenType generatorType = MMRG;
	ZZ m = power_ZZ(2,29)-3;
	int k = 3; 
	ZZ d = conv<ZZ>(0);
	ZZ c = conv<ZZ>(1);
	ZZ b; b = 2; b -= 2*c;
	BMat A;
	SavvidyMatrix(mat_ZZ_p& A, ZZ m, int k, ZZ d, ZZ c, ZZ b)
	//1
	int MinDim = 4;
	int MaxDim = 8;
	bool Dual = true;
	LatticeType latticeType = FULL;
	int Lacunary = 0;
	int MaxNodesBB = 1000000000;
	bool invertLength = false; // If true, the inverse of length of the shortest vector is printed, otherwise the length itself is printed.
	int details = 0;
	string ResultForm = "terminal";



	// building the MMRGLattice and LatTestSpectral objects



	// appliquer le test spectral


	

	// sortir les résultats



	return 0;
}
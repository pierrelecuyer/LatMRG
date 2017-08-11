#include <iostream>
#include "latmrg/LatConfig.h"
#include "latmrg/Writer.h"
#include "latmrg/WriterRes.h"
#include "latmrg/ReportLat.h"
#include "latmrg/ReportHeaderLat.h"
#include "latmrg/ReportFooterLat.h"
#include "latmrg/LatTestAll.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestBeyer.h"
#include "latmrg/MMRGLattice.h"

using namespace std;
using namespace LatticeTester;
using namespace LatMRG;

BMat SavvidyMatrix(BScal p, int N, BScal s, BScal m, BScal b) {
	ZZ_p::init(p);
	mat_ZZ_p A;
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

    return conv<BMat>(A);
}

//==========================================================================


//PW_TODO : clarifier notations

int main ()
{
	// parameters of the test
	//----------------------------------------------------------------------
	CriterionType criter = SPECTRAL;
	NormaType norma = BESTLAT;		
	NormType norm = L2NORM;					
	int J = 1;
	GenType genType = MMRG;
	//BScal m = power_ZZ(2,29)-3;
   BScal m = conv<ZZ>("2305843009213693951");
	int k = 8;
	BScal param_d = conv<ZZ>(0);
   BScal param_c = power_ZZ(2,53) + 1;
	BScal param_b; param_b = 2; param_b -= 2*param_c;
	BMat A = SavvidyMatrix(m, k, param_d, param_c, param_b);
	//1
	int d = 1;
	int fromDim = 9;
	int toDim = 38;
	bool dualF = true;
	LatticeType latType = FULL;

	int lacunary;
	int lacGroupSize = 0;
	if (lacGroupSize == 0)
		lacunary = 0;
	else
		lacunary = 1;
	BVect Lac;
	Lac.resize(lacGroupSize);
	//Lac[0] = 6;
	BScal lacSpacing;
	long maxNodesBB = 1000000000;
	bool invertF = false; // If true = inverse of length, if false = length itself
	int detailF = 0;
	OutputType outputType = TERMINAL;


	// creating a LatConfig object to store all the parameters
	//----------------------------------------------------------------------
	LatConfig config;

	config.criter = criter;
	if (config.criter == SPECTRAL)
  		config.norma = norma;
	
	config.norm = norm;
  	//config.readGenFile = readGenFile;
  	//if (config.readGenFile)
    	//readString (config.fileName, ln, 1);
    config.J = J;
	config.setJ (config.J);

	for (int i = 0; i < config.J; i++) {
		config.genType[i] = genType;
		config.comp[i] = new MRGComponent(m, A, k);
	}

	config.d = d;
	if (config.d < 1) {
	  MyExit (1, "ParamReaderLat:   config.d < 1");
	}
	config.td = new int[config.d];
	config.td[0] = fromDim;
	config.td[1] = toDim;


	config.dualF = dualF;
	config.latType = latType;

	if ((config.genType[0] == RANK1 || config.genType[0] == KOROBOV) && config.latType != FULL)
		MyExit (1, "ParamReaderLat:   latType must be FULL for KOROBOV or RANK1 lattices");
	 
	if (config.latType == ORBIT) {
		MyExit (1, "case ORBIT is not finished");
		//PW_TODO
	   	//readOrbit (config.J, config.comp, ++ln);
	}

	config.lacunary = lacunary;
	config.lacGroupSize = lacGroupSize;
	config.lacSpacing = lacSpacing;
	config.Lac = Lac;
	config.maxNodesBB = maxNodesBB;
	config.invertF = invertF;
	config.detailF = detailF;
	config.outputType = outputType;


	// building the MMRGLattice and LatTestSpectral objects
	//----------------------------------------------------------------------
   	string infile = "name";
    LatTestAll latTestAll;
    Writer* rw = latTestAll.createWriter (infile.c_str(), config.outputType);

	LatMRG::IntLattice *lattice = 0;
	LatMRG::IntLattice *master = 0;
	Lacunary *plac = 0;
	bool stationary = true;
	bool memLacF = true; 


	if (memLacF && config.lacunary) {
		lattice = new MMRGLattice (config.comp[0]->getM(), config.comp[0]->A,
									toDim,config.comp[0]->k, config.Lac, config.norm);           
	} else {
		lattice = new MMRGLattice (config.comp[0]->getM(), config.comp[0]->A,
									toDim,config.comp[0]->k, config.norm);
	}

	ReportHeaderLat header (rw, &config, lattice);
	ReportFooterLat footer (rw);
	ReportLat report (rw, &config, &header, &footer);

	double minVal[1 + toDim];
	SetZero (minVal, toDim);

	Normalizer *normal = 0;
	// creates and returns the normalizer corresponding to config.norma
	normal = lattice->getNormalizer (config.norma, 0, config.dualF); 
	normal->setNorm (config.norm);

	if (!memLacF && config.lacunary) {
		plac = new Lacunary (config.Lac, toDim);
		lattice->setLac (*plac);
	}

	switch (config.criter) {
		case SPECTRAL: {
				LatTestSpectral spectralTest (normal, lattice);
				lattice->buildBasis (fromDim - 1);
				spectralTest.attach (&report);
				report.printHeader ();
				spectralTest.setDualFlag (config.dualF);
				spectralTest.setInvertFlag (config.invertF);
				spectralTest.setDetailFlag (config.detailF);
				spectralTest.setMaxAllDimFlag (true);
				spectralTest.setMaxNodesBB (config.maxNodesBB);
				spectralTest.test (fromDim, toDim, minVal);
				// lattice->write();
				footer.setLatticeTest (&spectralTest);
				report.printTable ();
				report.printFooter ();
			}
			break;

		case BEYER: {
				LatTestBeyer beyerTest (lattice);
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

		default:
			cerr << "Default case for config.criter" << endl;
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


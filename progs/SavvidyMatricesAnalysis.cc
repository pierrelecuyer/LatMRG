#include <iostream>
#include "latmrg/LatTestSpectral.h"

using namespace std;
using namespace LatticeTester;
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
	//----------------------------------------------------------------------
	CriterionType criter = SPECTRAL;
	NormaType norma = BESTLAT;		
	NormType norm = L2NORM;
	//bool readGenFile = false;						
	int J = 1;

	MRGComponent **comp;

	GenType genType = MMRG; //GenType *genType;
	ZZ m = power_ZZ(2,29)-3;
	int k = 3; 
	ZZ d = conv<ZZ>(0);
	ZZ c = conv<ZZ>(1);
	ZZ b; b = 2; b -= 2*c;
	BMat A;
	SavvidyMatrix(mat_ZZ_p& A, ZZ m, int k, ZZ d, ZZ c, ZZ b)
	//1
	int d = 1;
	int *td;
	int fromDim = 4;
	int toDim = 8;
	bool dualF = true;
	LatticeType latType = FULL;

	int Lacunary = 0; // ou bool ?
	int lacGroupSize;
	BScal lacSpacing;
	BVect Lac;

	long MaxNodesBB = 1000000000;
	bool invertF = false; // If true, the inverse of length of the shortest vector is printed, otherwise the length itself is printed.
	int detailF = 0;
	OutputType outputType = "terminal";


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
	config.genType[i] = &genType;
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
	string fname (infile);
  	fname += ".dat";
   	ParamReaderLat paramRdr (fname.c_str ());
   	fname.clear ();
	Writer* rw = createWriter (infile, config.outputType);

	LatMRG::IntLattice *lattice = 0;
	LatMRG::IntLattice *master = 0;
	Lacunary *plac = 0;
	bool stationary = true;
	int toDim = config.td[1];
	int fromDim = config.td[0];
	bool memLacF = true; 


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
				if (1 == config.d) {
					spectralTest.test (fromDim, toDim, minVal);
					// lattice->write();
					footer.setLatticeTest (&spectralTest);
					report.printTable ();
					report.printFooter ();
				} else {
					if (config.genType[0] == MRG || config.genType[0] == LCG)
						master = new MRGLattice (*(MRGLattice *) lattice);
					else if (config.genType[0] == KOROBOV)
						master = new KorobovLattice (*(KorobovLattice *) lattice);
					else if (config.genType[0] == RANK1)
						master = new Rank1Lattice (*(Rank1Lattice *) lattice);

					master->buildBasis (toDim);
					TestProjections proj (master, lattice, &spectralTest, config.td, config.d);
					proj. setOutput (rw);
					// proj.setDualFlag (config.dualF);
					proj.setPrintF (true);
					double merit = proj.run (stationary, false, minVal);
					int nbProj = proj.getNumProjections ();
					rw->writeString ("\nMin merit:   ");
					rw->writeDouble (sqrt (merit));
					rw->newLine ();
					rw->writeString ("Num projections:   ");
					rw->writeInt (nbProj);
					rw->newLine ();
					// nbProj = proj.calcNumProjections(stationary, false);
					// cout << "Num projections2:  " << nbProj << endl << endl;
					delete master;
				}
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


	return 0;
}

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
  NTL::ZZ_p::init(p);
  NTL::mat_ZZ_p A;
  A.SetDims(N,N);
  for (int j = 1; j < N; j ++) {
    for (int i = j+1; i < N; i++)
      A[i][j] = NTL::conv<NTL::ZZ_p>((i-j+2) * m + b);
  }
  for (int i = 0; i < N; i ++)
    A[i][0] = 1;
  for (int i = 1; i < N; i++)
    A[i][i] = 2;
  for (int i = 0; i < N; i ++) {
    for (int j = i+1; j < N; j ++)
      A[i][j] = 1;
  }
  A[2][1] += NTL::conv<NTL::ZZ_p>(s);

  return NTL::conv<BMat>(A);
}

//==========================================================================



int main ()
{
  // parameters of the test
  //----------------------------------------------------------------------
  CriterionType criter = SPECTRAL;
  NormaType norma = BESTLAT;		
  NormType norm = L2NORM;					
  int J = 1;
  GenType genType = MMRG;
  BScal m = NTL::power_ZZ(2,61)-1;
  //BScal m = NTL::conv<NTL::ZZ>("2305843009213693951");
  int k = 240;
  BScal param_d = NTL::conv<NTL::ZZ>("487013230256099140");
  BScal param_c = NTL::power_ZZ(2,51) + 1;
  BScal param_b; param_b = 2; param_b -= 2*param_c;
  BMat A = SavvidyMatrix(m, k, param_d, param_c, param_b);

  //1
  int d = 1;
  int fromDim = 7;
  int toDim = 7;
  bool dualF = true;
  LatticeType latType = FULL;

  int lacunary;
  int lacGroupSize = 7;
  if (lacGroupSize == 0)
    lacunary = 0;
  else
    lacunary = 1;
  LacunaryType lactype(ARBITRARYINDICES);
  BVect Lac;
  Lac.resize(lacGroupSize);
  Lac[0] = 4;
  Lac[1] = 5;
  Lac[2] = 6;
  Lac[3] = 243;
  Lac[4] = 244;
  Lac[5] = 245;
  Lac[6] = 246;
  BScal lacSpacing;
  long maxNodesBB = 1000000000;
  bool invertF = false; // If true = inverse of length, if false = length itself
  int detailF = 0;
  OutputType outputType = TERMINAL;


  // creating a LatConfig object to store all the parameters
  //----------------------------------------------------------------------
  LatConfig<MScal> config;

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
    config.comp[i] = new MRGComponent<MScal>(m, A, k);
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
    //readOrbit (config.J, config.comp, ++ln);
  }

  config.lacunary = lacunary;
  config.lacunaryType = lactype;
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

  LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *lattice = 0;
  Lacunary<BScal, BVect> *plac = 0;
  bool memLacF = true; 


  if (memLacF && config.lacunary) {
    lattice = new MMRGLattice<MScal> (config.comp[0]->getM(), config.comp[0]->A,
        toDim, config.comp[0]->k, config.lacunaryType, config.Lac, config.norm);           
  } else {
    lattice = new MMRGLattice<MScal> (config.comp[0]->getM(), config.comp[0]->A,
        toDim, config.comp[0]->k, config.norm);
  }

  ReportHeaderLat header (rw, &config, lattice);
  ReportFooterLat footer (rw);
  ReportLat report (rw, &config, &header, &footer);

  double minVal[1 + toDim];
  SetZero (minVal, toDim);

  Normalizer<RScal> *normal = 0;
  // creates and returns the normalizer corresponding to config.norma
  normal = lattice->getNormalizer (config.norma, 0, config.dualF); 
  normal->setNorm (config.norm);

  if (!memLacF && config.lacunary) {
    plac = new Lacunary<BScal, BVect> (config.Lac, toDim);
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
                  LatTestBeyer<MScal, BScal, BVect, BMat, NScal, NVect, RScal> beyerTest (lattice);
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

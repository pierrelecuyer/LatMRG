/*
 * This file showcases the usage of the LatMRG API at a low level. It contains
 * the computations of the MIXMAX article (to be published). This program does 
 * something very similar to HighMixMax.cc but takes its parameters from the 
 * standard input instead of having hard coded parameters.
 *
 * To compile this program you should use something similar to :
 *      g++ -I./include -I./latticetester/include -c \
 *      ./examples/LowMixMax.cc -L./lib -llatmrg -llatticetester \
 *      -L/usr/local/lib -lntl -lgmp -o ./bin/examples/LowMixMax
 * with LatMRG as working directory with the LatMRG and LatticeTester library 
 * built.
 *
 * You can use this program with whatever parameters you see fit. It takes
 * either 6 or 7 arguments. The first one indicates how the MixMax generator 
 * will be constructed. It can be either 4 or 5, the number of parameters for
 * the construction. The next 4 or 5 parameters are the ones used to construct 
 * the generator (see the article for some precision on what each parameter
 * means. The last one indicates which projections to use. It can be 5, 7 or 0.
 * 5 tests the case for proposition 4 in the article, 7 tests proposition 5 in
 * the article and 0 tests indices in order from dimension 1 to dimension 48.
 * You can get the tests in the MixMax article by calling :
 *      ./bin/examples/LowMixMax < ./examples/MixMaxLatticeAnalysis/xx.txt
 * in the LatMRG directory (assuming you compiled with the command above, xx are
 * the numbers of one of the files).
 *
 * This program is not super well done because it reads directly from cin but it 
 * works just fine for the application it implements.
 *
 * Authors: - Paul Wambergue
 *          - Marc-Antoine Savard
 * */

// We need to handle large integers, but we don't need to be that precise when 
// we do floating point
#define NTL_TYPES_CODE 2

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
using namespace LatMRG;

/*
 * A small function that allows us to build the generating matrix of a MMRG
 * of the MixMax type as described in "Spectrum and Entropy of C-systems. MIXMAX 
 * Random Number Generator" (Savvidy, Savvidy 2016).
 * */
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

//=============================================================================

int main ()
{
  // Initialization of the MixMax generating matrix
  BMat A;
  BScal m;
  int k;
  string s;
  std::cin >> k;
  if (k == 4) {
    // This is the case with 4 parameters
    cin >> s;
    m = NTL::conv<NTL::ZZ>(s.c_str());
    cin >> k;
    cin >> s;
    BScal param_d = NTL::conv<NTL::ZZ>(s.c_str());
    cin >> s;
    BScal param_c = NTL::conv<NTL::ZZ>(s.c_str());
    BScal param_b; param_b = 2; param_b -= 2*param_c;
    A = SavvidyMatrix(m, k, param_d, param_c, param_b);
  } else if (k == 5) {
    // This is the case with 5 parameters
    cin >> s;
    m = NTL::conv<NTL::ZZ>(s.c_str());
    cin >> k;
    cin >> s;
    BScal param_d = NTL::conv<NTL::ZZ>(s.c_str());
    cin >> s;
    BScal param_c = NTL::conv<NTL::ZZ>(s.c_str());
    cin >> s;
    BScal param_b = NTL::conv<NTL::ZZ>(s.c_str());
    A = SavvidyMatrix(m, k, param_d, param_c, param_b);
  }

  //---------------------------------------------------------------------------
  // Properties needed to instanciate the generator

  // MixMax generators are MMRG
  GenType genType = MMRG;
  // The number of MMRG components. MixMax only has one.
  int J = 1;
  // The test we are going to do.
  LatticeTester::CriterionType criter = LatticeTester::SPECTRAL;
  // The normalization used to get the figure of merit.
  LatticeTester::NormaType norma = LatticeTester::BESTLAT;		
  // The norm we used for computations
  LatticeTester::NormType norm = LatticeTester::L2NORM;

  // The types of projections tested. d=1 is the classic case where the 
  // projections on dimensions from fromDim to toDim (inclusively) will be tested
  int d = 1;
  // I have to understand that. This causes Segmentation faults in the basis 
  // construction for whatever reason.
  int fromDim = 1;
  int toDim = 0;
  // We perform the spectral test on the dual, meaning we will find the shortest
  // vector of the dual l and compute a figure of merit based on 1/l.
  bool dualF = true;
  // The type of lattice we will work with. FULL means we will build the lattice
  // for all states the generator can get to, no matter what the seed is.
  LatticeType latType = FULL;

  // We build the sets for which we will do the projections.
  // I don't know how to build these exactly. I have to had functionnalities to
  // automatically build the sets of the propositions of the article.
  int lacunary = 0;
  int lacGroupSize;
  cin >> lacGroupSize;
  BVect Lac;
  Lac.resize(lacGroupSize);
  if (lacGroupSize == 0) {
    // We test all the sets we can with indices in order
    toDim = 48;
  } else if (lacGroupSize == 5) {
    // This is the case of proposition 4 in the mixmax article
    lacunary = 1;
    Lac[0] = 4;
    Lac[1] = 5;
    Lac[2] = k+3;
    Lac[3] = k+4;
    Lac[4] = k+5;
    toDim = 5;
  } else if (lacGroupSize == 7) {
    // This is the case of proposition 5 in the mixmax article
    lacunary = 1;
    Lac[0] = 4;
    Lac[1] = 5;
    Lac[2] = 6;
    Lac[3] = k+3;
    Lac[4] = k+4;
    Lac[5] = k+5;
    Lac[6] = k+6;
    toDim = 7;
  }
  LacunaryType lactype(ARBITRARYINDICES);


  BScal lacSpacing;
  // The maximum number of nodes in the BB algorithm to find the shortest vector
  long maxNodesBB = 1000000000;
  // If this is true, the inverse of the length of the shortest vector will be 
  // printed
  bool invertF = false;
  // The amount of detail printed. 0 is the least and the default.
  int detailF = 0;
  // We will print on the terminal.
  LatticeTester::OutputType outputType = LatticeTester::TERMINAL;

  //---------------------------------------------------------------------------
  // Creating a LatConfig object to store all the parameters
  // This section basicaly copies all the parameters we just defined in a 
  // LatConfig object.

  LatConfig<MScal> config;

  config.criter = criter;
  if (config.criter == LatticeTester::SPECTRAL)
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
    LatticeTester::MyExit (1, "ParamReaderLat:   config.d < 1");
  }
  config.td = new int[config.d];
  config.td[0] = fromDim;
  config.td[1] = toDim;


  config.dualF = dualF;
  config.latType = latType;

  if ((config.genType[0] == RANK1 || config.genType[0] == KOROBOV) && config.latType != FULL)
    LatticeTester::MyExit (1, "ParamReaderLat:   latType must be FULL for KOROBOV or RANK1 lattices");

  if (config.latType == ORBIT) {
    LatticeTester::MyExit (1, "case ORBIT is not finished");
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


  //---------------------------------------------------------------------------
  // building the MMRGLattice and LatTestSpectral objects

  string infile = "name";
  LatTestAll latTestAll;
  Writer* rw = latTestAll.createWriter (infile.c_str(), config.outputType);

  LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *lattice = 0;
  LatticeTester::Lacunary<BScal, BVect> *plac = 0;
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
  LatticeTester::SetZero (minVal, toDim);

  LatticeTester::Normalizer<RScal> *normal = 0;
  // creates and returns the normalizer corresponding to config.norma
  normal = lattice->getNormalizer (config.norma, 0, config.dualF); 
  normal->setNorm (config.norm);

  if (!memLacF && config.lacunary) {
    plac = new LatticeTester::Lacunary<BScal, BVect> (config.Lac, toDim);
    lattice->setLac (*plac);
  }

  switch (config.criter) {
    case LatticeTester::SPECTRAL: {
                                    LatTestSpectral spectralTest (normal, lattice);
                                    lattice->buildBasis (fromDim);
                                    spectralTest.attach (&report);
                                    //report.printHeader ();
                                    spectralTest.setDualFlag (config.dualF);
                                    spectralTest.setInvertFlag (config.invertF);
                                    spectralTest.setDetailFlag (config.detailF);
                                    spectralTest.setMaxAllDimFlag (true);
                                    spectralTest.setMaxNodesBB (config.maxNodesBB);
                                    spectralTest.test (fromDim, toDim, minVal);
                                    lattice->write();
                                    footer.setLatticeTest (&spectralTest);
                                    report.printTable ();
                                    report.printFooter ();
                                  }
                                  break;

    case LatticeTester::BEYER: {
                                 LatTestBeyer<MScal, BScal, BVect, BMat, NScal, NVect, RScal> beyerTest (lattice);
                                 lattice->buildBasis (fromDim);
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

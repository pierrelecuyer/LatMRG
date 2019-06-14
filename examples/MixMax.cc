#include "MixMaxConfigs.h"

#define LATMRG_USE_CONFIG
#define LATMRG_USE_TEST
#include "latmrg/Test.h"
#include "latmrg/MixmaxMMRG.h"
#include "latmrg/MMRGLattice.h"

using namespace LatMRG;

//void doTests(MScal m, int k, MScal d, MScal c) {
//  MixmaxMMRG<MScal> mixmax (m, k, d, c);
//  //cout << "A = \n" << mixmax.getMatrix() << endl;
//
//  LatConfig<MScal> latconfig;
//
//  latconfig.J = 1;
//  latconfig.setJ(latconfig.J);
//  for (int i = 0; i < latconfig.J; i++) {
//    latconfig.genType[i] = MMRG; //PW_TODO moche
//    latconfig.comp[i] = new MRGComponent<MScal>(m, mixmax.getMatrix(), k);
//  }
//
//  //latconfig.comp[0] = new MRGComponent (m, A, k);
//  //latconfig.genType[0] = MMRG;
//
//  latconfig.d = 1;
//  latconfig.td = new int[latconfig.d];
//
//  latconfig.criter = SPECTRAL;
//  latconfig.norma = BESTLAT;
//  latconfig.norm = L2NORM;
//  latconfig.dualF = true;
//  latconfig.invertF = false;
//  latconfig.latType = LatMRG::FULL;
//  latconfig.lacunary = true;
//  latconfig.lacunaryType = ARBITRARYINDICES;
//  latconfig.numberLacIndices = 7;
//  latconfig.maxNodesBB = 1000000;
//
//  // Proposition 1
//  latconfig.td[0] = 3;
//  latconfig.td[1] = 3;
//  latconfig.Lac.resize(3);
//  latconfig.Lac[0]=1;
//  latconfig.Lac[1]=k;
//  latconfig.Lac[2]=k+1;
//
//  std::vector<double> FoM;
//  FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);
//
//  printResult(FoM, latconfig.td[0]);
//
//  // Proposition 4
//  latconfig.td[0] = 5;
//  latconfig.td[1] = 5;
//  latconfig.Lac.resize(5);
//  latconfig.Lac[0]=4;
//  latconfig.Lac[1]=5;
//  latconfig.Lac[2]=k+3;
//  latconfig.Lac[3]=k+4;
//  latconfig.Lac[4]=k+5;
//
//  FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);
//
//  printResult(FoM, latconfig.td[0]);
//
//  // Proposition 6
//  latconfig.td[0] = 7;
//  latconfig.td[1] = 7;
//  latconfig.Lac.resize(7);
//  latconfig.Lac[0]=4;
//  latconfig.Lac[1]=5;
//  latconfig.Lac[2]=6;
//  latconfig.Lac[3]=k+3;
//  latconfig.Lac[4]=k+4;
//  latconfig.Lac[5]=k+5;
//  latconfig.Lac[6]=k+6;
//
//  FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);
//
//  printResult(FoM, latconfig.td[0]);
//
//  // Getting a bound on M_s
//  latconfig.d = 10;
//  delete[] latconfig.td;
//  latconfig.lacunary = false;
//  latconfig.td = new int[latconfig.d];
//  latconfig.td[0] = 10;
//  latconfig.td[1] = 10;
//  latconfig.td[2] = 10;
//  latconfig.td[3] = 10;
//  latconfig.td[4] = 10;
//  latconfig.td[5] = 10;
//  latconfig.td[6] = 10;
//  latconfig.td[7] = 10;
//  latconfig.td[8] = 10;
//  latconfig.td[9] = 10;
//
//  FoM = ComputeFigureOfMerit<MScal, BScal, BVect, BMat, NScal, NVect, RScal>(latconfig);
//
//  printResult(FoM, latconfig.td[0]);
//
//}

//==================================================================================

int main ()
{

  Config conf;
  conf.type = MMRG;
  conf.criterion = LatticeTester::SPECTRAL;
  conf.reduction = LatticeTester::FULL;
  conf.use_dual = true;
  for (int i = 0; i < num_confs; i++) {
    MixmaxMMRG<Int> matrix (Configs[i].m, Configs[i].k, Configs[i].d, Configs[i].c);
    conf.modulo = Configs[i].m;
    conf.order = Configs[i].k;
    {
      std::vector<std::size_t> max_coords;
      max_coords.push_back(0);
      max_coords.push_back(10);
      max_coords.push_back(Configs[i].k);
      Projections proj(3, 1, max_coords);
      conf.proj = &proj;
      conf.numProj = 3;
      conf.minDim = 1;
      conf.maxDim = 0;
      conf.projDim = max_coords;

      MMRGLattice<Int, Dbl> mixmax(Configs[i].m, matrix.getMatrix(), Configs[i].k+2, Configs[i].k);
      auto FoM = test(mixmax, conf);
      FoM.computeMerit("min");

      std::cout << "merit: " << FoM.getMerit() << "\n";
      std::cout << "Worst proj: " << FoM.worstProj() << "\n";
      std::cout << "Short vector: " << FoM.worstVect() << "\n";
    }

  }

  return 0;
}

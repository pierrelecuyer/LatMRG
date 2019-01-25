/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "latticetester/NormaBestLat.h"
#include "latticetester/Const.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/ParamReaderExt.h"

using namespace LatMRG;

namespace {
  typedef NTL::ZZ Int;
  typedef NTL::RR Dbl;

  Chrono timer;
  long num_gen = 0;

  Int b = NTL::power2_ZZ(64);
  int minDim = 1, maxDim = 20;

  /**
   * Prints the results of the test on lattice stored in latTest.
   * */
  void printTest(MRGLattice<Int, Dbl>& lattice, Dbl minVal[]) {
    std::cout << lattice.getCoef() << std::endl;
    for (int i = minDim; i <= maxDim; i ++) {
      std::cout << "dim[" << i << "], shortest: " << minVal[i-minDim] << std::endl;
    }

  }

  /**
   * Reads the configuration file for this program.
   * */
  double test(MRGLattice<Int, Dbl> & lattice) {
    Dbl minVal[maxDim - minDim + 1];
    LatticeTester::Normalizer<Dbl>* norma(lattice.getNormalizer(LatticeTester::BESTLAT, 0, true));
    for (int i = minDim; i <= maxDim; i ++) {
      lattice.buildBasis(i);
      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lattice);
      lattice.dualize();
      red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lattice.getDim());
      red.shortestVector(lattice.getNorm());
      NTL::vector<Int> shortest(lattice.getBasis()[0]);
      Dbl tmp;
      LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
      minVal[i-minDim] = Dbl(1)/(NTL::sqrt(tmp));
      std::cout << "Normalized value: " << NTL::sqrt(tmp)/norma->getBound(i) << std::endl;
      lattice.dualize();
    }
    printTest(lattice, minVal);
    delete norma;
    return 0;
    //latTest.getMerit().getST(minDim, maxDim);
  }


  // This just instanciates number MWC generators with order k and mod b.
  void testGenerators() {
    b = 6;
    NTL::vector<Int> e;
    e.resize(22);
    for (int i = 0; i < 22; i++) e[i] = 0;
    e[2] = 1;
    e[21] = 1;
    MWCLattice<Int, Dbl> lattice(b, e, 21);
    Int a(lattice.getCoef()[0]);
    Int m(lattice.getModulo());
    NTL::vector<Int> a_vec(2);
    for (int i = 1; i < 20; i+=2) {
      std::cout << "i: " << i << std::endl;
      a_vec[1] = NTL::power(a, i)%m;
      MRGLattice<Int, Dbl> lattice2(m, a_vec, maxDim, 1, FULL);
      test(lattice2);
    }
    num_gen++;
  }
}

int main (int argc, char **argv)
{
  timer.init();
  testGenerators();
  std::cout << "CPU time: " << timer.toString() << std::endl;
  return 0;
}

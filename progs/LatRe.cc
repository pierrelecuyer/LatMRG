#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "latticetester/Types.h"

#include "latmrg/LatTestAll.h"

using namespace LatMRG;
using LatticeTester::IntLattice;

namespace {
  typedef NTL::ZZ Int;
  typedef double Dbl;
  typedef NTL::matrix<Int> IntMat;
  typedef NTL::vector<Int> IntVec;

  // Program global objects
  Chrono timer;
  LatticeTester::Normalizer<Dbl>* norma; // If a normalizer is used
  // For period tests
  DecompType decompm1 = DECOMP, decompr = DECOMP;
  std::string filem1, filer;
  Dbl worstMerit(1);

  // Data file read parameters
  GenType type;
  int minDim, maxDim;
  LatticeTester::NormaType normaType;
  double timeLimit;

  // MRG Specific parameters
  IntVec mult; // MRG multipliers

  // MWC Specific parameters
  Int b; // modulo of MWC recurence

  // Shared components names
  std::int64_t order; // k for MRG and MMRG
  Int modulo; // m for MRG and MMRG
  // modulo is basis^exponent+rest
  std::int64_t basis;
  std::int64_t exponent;
  std::int64_t rest;
  bool period;

  /**
   * Prints the results of the program execution.
   * */
  void printResults() {
    std::cout << "Merit: " << worstMerit << std::endl;
    std::cout << "CPU time: " << timer.toString() << std::endl;
  }

  /**
   * Tests the generator via spectral test.
   * */
  Dbl test(IntLattice<Int, Int, Dbl, Dbl> & lattice) {
    norma = lattice.getNormalizer(normaType, 0, true);
    Dbl merit = Dbl(1);
    for (int i = minDim; i <= maxDim; i++){
      //std::cout << "i: " << i << std::endl;
      // Building the basis
      lattice.buildBasis(i);
      lattice.dualize();
      // Reducing the lattice
      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lattice);
      red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lattice.getDim());
      //red.redLLLNTL(0.99, LatticeTester::QUADRUPLE, lattice.getDim());
      red.shortestVector(lattice.getNorm());
      // Computing shortest vector length and spectral test
      NTL::vector<Int> shortest(lattice.getBasis()[0]);
      Dbl tmp;
      LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
      // Normalization
      tmp = NTL::sqrt(tmp)/norma->getBound(i);
      // The next condition is true if there is no normalizer (or a bug)
      if (tmp > 1) tmp = Dbl(1)/tmp;
      //std::cout << "Value: " << tmp << std::endl;
      merit = (tmp < merit) ? tmp : merit;
    }

    delete norma;
    return merit;
  }

  // This is the main program function. This instanciates evey generator to
  // test and tests it.
  void testGenerator() {
    if (type == MRG) {
      IntVec temp(order+1);
      temp[0] = Int(0);
      for (int i = 1; i < order+1; i++) temp[i] = mult[i-1];
      MRGLattice<Int, Dbl> mrglat(modulo, temp, maxDim, order, FULL);
      worstMerit = test(mrglat);
    } else if (type == MWC) {
      //MWCLattice<Int, Dbl>* mwclat = 0;
      //mwclat = nextGenerator(mwclat);
      //Dbl merit(test(*mwclat));
      //if (merit > bestMerit) {
      //  addLattice(mwclat, merit);
      //}
    } else if (type == MMRG) {
      //MMRGLattice<Int, Dbl>* mmrglat = 0;
      // mmrglat = nextGenerator(mmrglat);
      // Dbl merit(test(*mmrglat));
      // if (merit > bestMerit) {
      //   addLattice(mmrglat, merit);
      // }
    }
  }

  // Reads parameters from config file.
  bool readConfigFile(int argc, char **argv) {
    ParamReaderExt<Int, Dbl> reader(argv[1]);
    reader.getLines();
    //int power;
    unsigned int ln = 0;
    reader.readGenType(type, ln++, 0);
    reader.readInt(minDim, ln, 0);
    reader.readInt(maxDim, ln++, 1);
    reader.readNormaType(normaType, ln++, 0);
    reader.readDouble(timeLimit, ln++, 0);
    if (type == MRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
      minDim = minDim<order ? order : minDim;
      maxDim = maxDim>minDim ? maxDim : minDim;
      mult.SetLength(order);
      reader.readMVect(mult, ln, 0, order, 0);
      ln++;
      reader.readBool(period, ln, 0);
      if (period) {
        // Using default parameters
        bool def;
        reader.readBool(def, ln, 1);
        if (!def) {
          reader.readDecompType(decompm1, ln, 2);
          reader.readString(filem1, ln, 3);
          reader.readDecompType(decompr, ln, 4);
          reader.readString(filer, ln, 5);
        }
      }
    } else if (type == MWC) {
      //reader.readInt(power, ln++, 0);
      //b = NTL::power2_ZZ(power);
      //reader.readLong(exponent, ln++, 0);
    } else if (type == MMRG) {
      //reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      //reader.readLong(order, ln++, 0);
    }
    return true;
  }

}

//==========================================================================

int main (int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return -1;
  }
  // Initializing values
  srand(time(NULL));
  filem1 = "./tempm1" + std::to_string(rand());
  filer = "./tempr" + std::to_string(rand());
  readConfigFile(argc, argv);
  timer.init();
  // Testing the generator(s)
  testGenerator();
  printResults();
  //delete bestLattice;
  //delete mod;
  return 0;
}


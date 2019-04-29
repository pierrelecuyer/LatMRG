#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "latticetester/Types.h"

#include "Exec.h"

using namespace LatMRG;
using LatticeTester::IntLattice;

namespace {
  // Program global objects
  Chrono timer;
  LatticeTester::Normalizer<Dbl>* norma; // If a normalizer is used
  // For period tests
  DecompType decompm1 = DECOMP, decompr = DECOMP;
  std::string filem1, filer;
  TestList* bestLattice;

  std::string delim = "\n========================================"
    "========================================\n\n";

  // Data file read parameters
  GenType type;
  int numProj;
  int minDim, maxDim;
  std::vector<std::size_t> projDim;
  Projections* proj;

  LatticeTester::NormaType normaType;
  double timeLimit;
  int num_gen;

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
    std::cout << "LatRe: A program to test Random Number Generators\n";
    std::cout << delim;
    std::cout << "Bellow are the results of the tests of " << num_gen << "generators:\n";
    std::cout << "Generator type: " << toStringGen(type) << "\n";
    if (type == MRG) {
      std::cout << "With modulus:   m = " << modulo << " = " << basis << "^"
        << exponent << (rest>0?"+":"") << (rest!=0?std::to_string(rest):"") << "\n";
      std::cout << "Of order:       k = " << order << "\n";
    } else if (type == MWC) {
    } else if (type == MMRG) {
    }
    std::cout << "And " << (period?"full":"any") << " period length\n";
    std::cout << "The test was:\n";
    std::cout << "Minimal " << (normaType==LatticeTester::NONE?"":"normalized ")
      << "shortest" << " non-zero vector length\n";
    if (normaType != LatticeTester::NONE) {
      std::cout << "Normalizer used: "
        << LatticeTester::toStringNorma(normaType) << "\n";
    }
    std::cout << "On dimensions and projections:\n";
    //proj->resetDim();
    //while(!proj->end()) {
    //  std::cout << proj->next() << "\n";
    //}
    std::cout << delim;
    std::cout << "Allowed running time: " << timeLimit << "s.\n";
    std::cout << "Actual CPU time: " << timer.toString() << "\n";
    std::cout << "Number of generators tested: " << num_gen << "\n\n";
    std::cout << "Retained generators (from best to worst):\n";
    auto list = bestLattice->getList();
    for (auto it = list.begin(); it!= list.end(); it++) {
      std::cout << delim;
      if (type == MRG) {
        std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
      } else if (type == MWC) {}
      else if (type == MMRG) {}
      std::cout << "Merit: " << (*it).getMerit() << "\n";
    }
  }

  void reduce(IntLattice<Int, Int, Dbl, Dbl>& lattice) {
    LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lattice);
    red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lattice.getDim());
    red.shortestVector(lattice.getNorm());
  }

  /**
   * Tests the generator via spectral test.
   * */
  DblVec test(IntLattice<Int, Int, Dbl, Dbl> & lattice) {
    norma = lattice.getNormalizer(normaType, 0, true);
    DblVec results(0);
    lattice.buildBasis(minDim);
    for (int i = minDim; i <= maxDim; i++){
      // Building the basis
      lattice.dualize();
      // Reducing the lattice
      reduce(lattice);
      // Computing shortest vector length and spectral test
      NTL::vector<Int> shortest(lattice.getBasis()[0]);
      Dbl tmp;
      LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
      // Normalization
      tmp = NTL::sqrt(tmp)/norma->getBound(i);
      // The next condition is true if there is no normalizer (or a bug)
      if (tmp > 1) tmp = Dbl(1)/tmp;
      results.append(tmp);
      lattice.dualize();
      lattice.incDim();
    }

    // Testing projections if there are any
    for (int i = 2; i <= numProj; i++) {
      proj->resetDim(i);
      lattice.buildBasis(projDim[i-1]+1);
      while(!proj->end(1)) {
        // Building the projection
        IntLattice<Int, Int, Dbl, Dbl> proj_lat(modulo, order, i, true);
        LatticeTester::Coordinates iter(proj->next());
        lattice.buildProjection(&proj_lat, iter);
        norma->setLogDensity(Dbl(-i*log(modulo)
              +log(abs(NTL::determinant(proj_lat.getBasis())))));
        proj_lat.dualize();
        // Reduction
        reduce(proj_lat);
        // Figure of merit
        IntVec shortest(proj_lat.getBasis()[0]);
        Dbl tmp;
        LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
        tmp = NTL::sqrt(tmp)/norma->getBound(i);
        if (tmp > 1) tmp = Dbl(1)/tmp;
        results.append(tmp);
      }
    }

    delete norma;
    return results;
  }

  // This is the main program function. This instanciates evey generator to
  // test and tests it.
  void testGenerator() {
    if (type == MRG) {
      IntVec temp(order+1);
      temp[0] = Int(0);
      for (int i = 1; i < order+1; i++) temp[i] = mult[i-1];
      MRGLattice<Int, Dbl> mrglat(modulo, temp, maxDim, order, FULL);
      Test the_test(mrglat.toString(), test(mrglat));
      bestLattice->add(the_test);
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
    reader.readInt(numProj, ln++, 0);
    reader.readInt(minDim, ln, 0);
    reader.readInt(maxDim, ln++, 1);
    projDim.push_back((unsigned)(maxDim-1));
    for (int i = 0; i<numProj-1; i++) {
      int tmp;
      reader.readInt(tmp, ln, i);
      // If the projection is requested only on indices smaller or equal to the
      // dimension, we learn nothing, this is corrected
      tmp = (unsigned)tmp>projDim.size()+1?tmp:projDim.size()+1;
      projDim.push_back((unsigned)(tmp-1));
    }
    ln++;
    reader.readNormaType(normaType, ln++, 0);
    reader.readInt(num_gen, ln++, 0);
    reader.readDouble(timeLimit, ln++, 0);
    if (type == MRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
      // Making sure that minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      mult.SetLength(order);
      reader.readMVect(mult, ln, 0, order, 0);
      ln++;
      reader.readBool(period, ln, 0);
      minDim = (!period||minDim>order) ? minDim : (order+1);
      maxDim = maxDim>minDim ? maxDim : minDim;
      for (unsigned int i = 0; i < projDim.size(); i++) {
        projDim[i] = (projDim[i]<(unsigned)(minDim-1))?(unsigned)(minDim-1):projDim[i];
      }
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
  bestLattice = new TestList(num_gen);
  proj = new Projections(numProj, minDim, projDim);
  timer.init();
  // Testing the generator(s)
  testGenerator();
  printResults();
  delete bestLattice;
  delete proj;
  return 0;
}


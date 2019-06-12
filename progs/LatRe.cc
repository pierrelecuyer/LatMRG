/**
 * This is a program that tests random number generators specified in a file
 *
 * Todo
 * - Test period
 * - Read from generator file
 * - input file format correction
 * */

#define LATMRG_TEST
#include "Exec.h"

using namespace LatMRG;
using namespace LatticeTester::Random;

namespace {
  Config conf;
  // Program global objects
  Chrono timer;

  /**
   * Prints the results of the program execution.
   * */
  void printResults() {
    std::cout << "LatRe: A program to test Random Number Generators\n";
    std::cout << delim;
    std::cout << "Bellow are the results of the tests of " << conf.num_gen << " generators:\n";
    std::cout << "Generator type: " << toStringGen(conf.type) << "\n";
    if (conf.type == MRG) {
      std::cout << "With modulus:   m = " << conf.modulo << " = " << conf.basis << "^"
        << conf.exponent << (conf.rest>0?"+":"") << (conf.rest!=0?std::to_string(conf.rest):"") << "\n";
      std::cout << "Of order:       k = " << conf.order << "\n";
    } else if (conf.type == MWC) {
    } else if (conf.type == MMRG) {
    }
    std::cout << "And " << (conf.period?"full":"any") << " period length\n";
    std::cout << "The test was:\n";
    std::cout << "Minimal " << (conf.normaType==LatticeTester::NONE?"":"normalized ")
      << "shortest" << " non-zero vector length\n";
    if (conf.normaType != LatticeTester::NONE) {
      std::cout << "Normalizer used: "
        << LatticeTester::toStringNorma(conf.normaType) << "\n";
    }
    std::cout << "On dimensions and projections:\n";
    std::cout << conf.proj->toString();
    std::cout << delim;
    std::cout << "Allowed running time: " << conf.timeLimit << "s.\n";
    std::cout << "Actual CPU time: " << timer.toString() << "\n";
    std::cout << "Number of generators tested: " << conf.num_gen << "\n\n";
    auto list = conf.bestLattice->getList();
    for (auto it = list.begin(); it!= list.end(); it++) {
      std::cout << delim;
      if (conf.type == MRG) {
        std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
      } else if (conf.type == MWC) {}
      else if (conf.type == MMRG) {}
      std::cout << "Merit: " << (*it).getMerit() << "\n";
    }
  }


  // This is the main program function. This instanciates evey generator to
  // test and tests it.
  void testGenerator() {
    if (conf.type == MRG) {
      IntVec temp(conf.order+1);
      temp[0] = Int(0);
      for (int i = 1; i < conf.order+1; i++) temp[i] = conf.mult[i-1];
      MRGLattice<Int, Dbl> mrglat(conf.modulo, temp, conf.maxDim, conf.order, FULL);
      Test the_test(mrglat.toString(), test(mrglat, conf));
      conf.bestLattice->add(the_test);
    } else if (conf.type == MWC) {
      //MWCLattice<Int, Dbl>* mwclat = 0;
      //mwclat = nextGenerator(mwclat);
      //Dbl merit(test(*mwclat));
      //if (merit > bestMerit) {
      //  addLattice(mwclat, merit);
      //}
    } else if (conf.type == MMRG) {
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
    reader.readGenType(conf.type, ln++, 0);
    reader.readCriterionType(conf.criterion, ln, 0);
    if (conf.criterion == LatticeTester::SPECTRAL) {
      reader.readBool(conf.use_dual, ln, 1);
      reader.readPreRed(conf.reduction, ln, 2);
      reader.readNormaType(conf.normaType, ln++, 3);
    } else if (conf.criterion == LatticeTester::LENGTH) {
      reader.readBool(conf.use_dual, ln, 1);
      reader.readPreRed(conf.reduction, ln++, 2);
    } else if (conf.criterion == LatticeTester::BEYER) {
      conf.reduction = LatticeTester::NOPRERED;
      ln++;
    }
    reader.readInt(conf.numProj, ln++, 0);
    reader.readInt(conf.minDim, ln, 0);
    reader.readInt(conf.maxDim, ln++, 1);
    conf.projDim.push_back((unsigned)(conf.maxDim-1));
    for (int i = 0; i<conf.numProj-1; i++) {
      int tmp;
      reader.readInt(tmp, ln, i);
      // If the projection is requested only on indices smaller or equal to the
      // dimension, we learn nothing, this is corrected
      tmp = (unsigned)tmp>conf.projDim.size()+1?tmp:conf.projDim.size()+1;
      conf.projDim.push_back((unsigned)(tmp-1));
    }
    ln++;
    reader.readInt(conf.num_gen, ln++, 0);
    reader.readDouble(conf.timeLimit, ln++, 0);
    if (conf.type == MRG) {
      reader.readNumber3(conf.modulo, conf.basis, conf.exponent, conf.rest, ln++, 0);
      reader.readLong(conf.order, ln++, 0);
      // Making sure that conf.minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      conf.mult.SetLength(conf.order);
      reader.readMVect(conf.mult, ln, 0, conf.order, 0);
      ln++;
      reader.readBool(conf.period, ln, 0);
      conf.minDim = (!conf.period||conf.minDim>conf.order) ? conf.minDim : (conf.order+1);
      conf.maxDim = conf.maxDim>conf.minDim ? conf.maxDim : conf.minDim;
      for (unsigned int i = 0; i < conf.projDim.size(); i++) {
        conf.projDim[i] = (conf.projDim[i]<(unsigned)(conf.minDim-1))?(unsigned)(conf.minDim-1):conf.projDim[i];
      }
      if (conf.period) {
        // Using default parameters
        bool def;
        reader.readBool(def, ln, 1);
        if (!def) {
          reader.readDecompType(conf.decompm1, ln, 2);
          reader.readString(conf.filem1, ln, 3);
          reader.readDecompType(conf.decompr, ln, 4);
          reader.readString(conf.filer, ln, 5);
        }
      }
    } else if (conf.type == MWC) {
      //reader.readInt(power, ln++, 0);
      //b = NTL::power2_ZZ(power);
      //reader.readLong(conf.exponent, ln++, 0);
    } else if (conf.type == MMRG) {
      //reader.readNumber3(modulo, conf.basis, conf.exponent, conf.rest, ln++, 0);
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
  conf.filem1 = "./tempm1" + std::to_string(rand());
  conf.filer = "./tempr" + std::to_string(rand());
  readConfigFile(argc, argv);
  conf.bestLattice = new TestList(conf.num_gen, true);
  conf.proj = new Projections(conf.numProj, conf.minDim, conf.projDim);
  timer.init();
  // Testing the generator(s)
  testGenerator();
  printResults();
  delete conf.bestLattice;
  delete conf.proj;
  return 0;
}

/**
 * This is a program that tests random number generators specified in a file
 *
 * Todo
 * - Test period
 * - Read from generator file
 * - input file format correction
 * */
#include "Exec.h"

using namespace LatMRG;
using namespace LatticeTester::Random;

namespace {
  Config conf;
  // Program global objects
  Chrono timer;
  int numProj, minDim, maxDim;
  int detail;
  std::vector<int> projDim;

  /**
   * Prints the results of the program execution.
   * */
  template<typename Lat>
    void printResults(MeritList<Lat>& bestLattice) {
      std::cout << "LatRe: A program to test Random Number Generators\n";
      std::cout << delim;
      std::cout << "Bellow are the results of the tests of " << conf.max_gen << " generators:\n";
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
      std::cout << "Number of generators tested: " << conf.max_gen << "\n\n";
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        std::cout << delim;
        if (conf.type == MRG) {
          std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
        } else if (conf.type == MWC) {}
        else if (conf.type == MMRG) {}
        if (detail == 0) {
          std::cout << (*it).toStringMerit();
        } else if (detail == 1) {
          std::cout << (*it).toStringProjections();
        }
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
    reader.readInt(numProj, ln++, 0);
    reader.readInt(minDim, ln, 0);
    reader.readInt(maxDim, ln++, 1);
    projDim.push_back((maxDim-1));
    for (int i = 0; i<numProj-1; i++) {
      int tmp;
      reader.readInt(tmp, ln, i);
      // If the projection is requested only on indices smaller or equal to the
      // dimension, we learn nothing, this is corrected
      if (tmp > 0 && tmp < i+3) {
        tmp = i+3;
      }
      projDim.push_back(tmp-1);
    }
    ln++;
    reader.readInt(conf.max_gen, ln++, 0);
    reader.readDouble(conf.timeLimit, ln++, 0);
    reader.readInt(detail, ln++, 0);
    if (conf.type == MRG) {
      reader.readNumber3(conf.modulo, conf.basis, conf.exponent, conf.rest, ln++, 0);
      reader.readLong(conf.order, ln++, 0);
      // Making sure that minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      conf.mult.SetLength(conf.order);
      reader.readMVect(conf.mult, ln, 0, conf.order, 0);
      ln++;
      reader.readBool(conf.period, ln, 0);
      minDim = (!conf.period||minDim>conf.order) ? minDim : (conf.order+1);
      maxDim = maxDim>minDim ? maxDim : minDim;
      for (unsigned int i = 0; i < projDim.size(); i++) {
        if (projDim[i]>0) projDim[i] = (projDim[i]<(minDim-1))?(minDim-1):projDim[i];
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
      reader.readNumber3(conf.modulo, conf.basis, conf.exponent, conf.rest, ln++, 0);
      reader.readLong(conf.order, ln++, 0);
      // Making sure that minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      conf.matrix.resize(conf.order, conf.order);
      reader.readBMat(conf.matrix, ln, 0, conf.order);
      reader.readBool(conf.period, ln, 0);
      minDim = (!conf.period||minDim>conf.order) ? minDim : (conf.order+1);
      maxDim = maxDim>minDim ? maxDim : minDim;
      for (unsigned int i = 0; i < projDim.size(); i++) {
        if (projDim[i] > 0) projDim[i] = (projDim[i]<(minDim-1))?(minDim-1):projDim[i];
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
  conf.proj = new Projections(numProj, minDim, projDim);
  timer.init();

  // Testing the generator(s)
  if (conf.type == MRG) {
    MeritList<MRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
    IntVec temp(conf.order+1);
    temp[0] = Int(0);
    for (int i = 1; i < conf.order+1; i++) temp[i] = conf.mult[i-1];
    MRGLattice<Int, Dbl> mrglat(conf.modulo, temp, maxDim, conf.order, FULL);
    bestLattice.add(test(mrglat, conf));
    printResults(bestLattice);
  } else if (conf.type == MWC) {
    //MWCLattice<Int, Dbl>* mwclat = 0;
    //mwclat = nextGenerator(mwclat);
    //Dbl merit(test(*mwclat));
    //if (merit > bestMerit) {
    //  addLattice(mwclat, merit);
    //}
  } else if (conf.type == MMRG) {
    MeritList<MMRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
    MMRGLattice<Int, Dbl> mmrglat(conf.modulo, conf.matrix, maxDim, conf.order);
    bestLattice.add(test(mmrglat, conf));
    printResults(bestLattice);
  }
  delete conf.proj;
  return 0;
}

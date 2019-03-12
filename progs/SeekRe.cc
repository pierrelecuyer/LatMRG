/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "latticetester/NormaBestLat.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/ParamReaderExt.h"

using namespace LatMRG;
using LatticeTester::IntLattice;

namespace {
  typedef NTL::ZZ Int;
  typedef NTL::RR Dbl;
  typedef NTL::matrix<Int> IntMat;

  // Program global objects
  Chrono timer;
  LatticeTester::Normalizer<Dbl>* norma;
  IntLattice<Int, Int, Dbl, Dbl>* bestLattice = NULL;
  Dbl bestMerit = Dbl(0);
  long num_gen = 0;

  // Data file parameters
  GenType type; // This one is experimental

  // MWC Specific parameters
  Int b; // modulo of MWC recurence

  // Shared components names
  std::int64_t order; // k for MRG and MMRG
  Int modulo; // m for MRG and MMRG
  // modulo is basis^exponent+rest
  std::int64_t basis;
  std::int64_t exponent;
  std::int64_t rest;

  // Non type dependent parameters
  int minDim, maxDim;
  double timeLimit;

  /**
   * A small class to search for modulus for MWC generators.
   * */
  class Modulus {
    public:

      /**
       * Search will begin at `2^e±c` and will increase of decrease depending on
       * `increase`.
       * */
      Modulus(long e, Int c, bool increase) {
        m = (Int(1)<<e) - 1;
        this->increase = increase;
        if (increase) m+=c;
        else m-=c;
      }

      /**
       * This function finds the next value for `m`. This can return 1 if
       * `increase` is `false` and `m` gets to 1.
       * */
      Int next() {
        while (m > 1) {
          nextM();
          LatticeTester::PrimeType status = LatticeTester::IntFactor<Int>::isPrime (m, KTRIALS);
          if (status == LatticeTester::PRIME || status == LatticeTester::PROB_PRIME) {
            if (1 == m % 4) continue;
            Int m1 = (m - 1)/2;
            status = LatticeTester::IntFactor<Int>::isPrime (m1, KTRIALS);
            if (status != LatticeTester::PRIME && status != LatticeTester::PROB_PRIME) continue;
            // For MWC we check 2^{(m-1)/2} \neq 1 mod m.
            NTL::ZZ_p::init(m);
            NTL::ZZ_p a = NTL::power(NTL::ZZ_p(2), m1);
            if (NTL::IsOne(a)) {
              continue;
            }
            return m;
          }
        }
        return Int(-1);
      }

    private:

      static const long KTRIALS = 200;

      /**
       * Last valid modulus found (or 2^e).
       * */
      Int m;

      /**
       * As passed to constructor.
       * */
      bool increase;

      /**
       * Increment/decrement m
       * */
      void nextM () {
        if (increase) m += 2;
        else m -= 2;
        if (0 == m % 5){
          if (increase) m += 2;
          else m -= 2;
        }
      }
  } /**mod*/;

  /*
   * The goal is to create this overload and to use it to switch generators
   * without requiring the use of switch statements.
   * */
  MRGLattice<Int, Dbl>* nextGenerator(MRGLattice<Int, Dbl>* lattice) {
    return 0;
  }

  MWCLattice<Int, Dbl>* nextGenerator(MWCLattice<Int, Dbl>* lattice) {
    Int m(0);
    long exp = exponent-1;
    // 63 bits at a time because NTL converts from SIGNED long
    while(exp > 0) {
      if (exp < 63) {
        m << exp;
        m += LatticeTester::RandBits(exp);
        exp -= exp;
      }
      m << 63;
      m += LatticeTester::RandBits(63);
      exp -= 63;
    }
    if ((m&1) == 1) m+=1;
    Modulus mod(exponent, m, true);
    if (lattice) delete lattice;
    return new MWCLattice<Int, Dbl>(b, mod.next());
  }

  /*
   * Goin' full random for now
   * */
  MMRGLattice<Int, Dbl>* nextGenerator(MMRGLattice<Int, Dbl>* lattice) {
    IntMat A;
    A.SetDims(order, order);
    NTL::clear(A);
    while (NTL::determinant(A) == 0) {
      for (long i = 0; i<order; i++) {
        for (long j = 0; j<order; j++) {
          Int a(0);
          long exp = exponent-1;
          // 63 bits at a time because NTL converts from SIGNED long
          while(exp > 0) {
            if (exp < 63) {
              a << exp;
              a += LatticeTester::RandBits(exp);
              exp -= exp;
            }
            a << 63;
            a += LatticeTester::RandBits(63);
            exp -= 63;
          }
          A[i][j] = a;
        }
      }
    }
    if (lattice) delete lattice;
    return new MMRGLattice<Int, Dbl>(modulo, A, maxDim, order);
  }

  /*
   * These next function add the tested lattices to the list of the best ones.
   * This only add the lattices that are good enough.
   * */
  void addLattice(MRGLattice<Int, Dbl>* lattice, Dbl& merit) {
    bestMerit = merit;
    if(bestLattice) delete bestLattice;
    bestLattice = new IntLattice<Int, Int, Dbl, Dbl>(*lattice);
  }
  void addLattice(MWCLattice<Int, Dbl>* lattice, Dbl& merit) {
    bestMerit = merit;
    if(bestLattice) delete bestLattice;
    bestLattice = new IntLattice<Int, Int, Dbl, Dbl>(*lattice);
  }
  void addLattice(MMRGLattice<Int, Dbl>* lattice, Dbl& merit) {
    bestMerit = merit;
    if(bestLattice) delete bestLattice;
    bestLattice = new IntLattice<Int, Int, Dbl, Dbl>(*lattice);
  }

  /*
   * These next function add the tested lattices to the list of the best ones.
   * This only add the lattices that are good enough.
   * */
  void printResults() {
    std::cout << "Number of generators tested: " << num_gen << std::endl;
    std::cout << "Best merit: " << bestMerit << std::endl;
    //std::cout << "Best lattice LCG coefficient: " << bestLattice->getCoef() << std::endl;
    std::cout << "Best lattice modulus: " << bestLattice->getModulo() << std::endl;
    std::cout << "CPU time: " << timer.toString() << std::endl;
  }

  /**
   * Tests the generator via spectral test.
   * */
  Dbl test(IntLattice<Int, Int, Dbl, Dbl> & lattice) {
    norma = lattice.getNormalizer(LatticeTester::BESTLAT, 0, true);
    Dbl merit = Dbl(1);
    for (int i = minDim; i <= maxDim; i++){
      std::cout << "i: " << i << std::endl;
      // Building the basis
      lattice.buildBasis(i);
      lattice.dualize();
      // Reducing the lattice
      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lattice);
      //red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lattice.getDim());
      red.redLLLNTL(0.99, LatticeTester::QUADRUPLE, lattice.getDim());
      //red.shortestVector(lattice.getNorm());
      // Computing shortest vector length and spectral test
      NTL::vector<Int> shortest(lattice.getBasis()[0]);
      Dbl tmp;
      LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
      //std::cout << "Pre-Normalization: " << NTL::sqrt(tmp) << std::endl;
      // Normalization
      //std::cout << "Bound: " << norma->getBound(i) << std::endl;
      tmp = NTL::sqrt(tmp)/norma->getBound(i);
      //tmp = Dbl(1)/NTL::sqrt(tmp);
      std::cout << "Value: " << tmp << std::endl;
      merit = (tmp < merit) ? tmp : merit;
    }

    delete norma;
    return merit;
  }


  // This just instanciates number MWC generators with order k and mod b.
  void testGenerators() {
    MRGLattice<Int, Dbl>* mrglat = 0;
    MMRGLattice<Int, Dbl>* mmrglat = 0;
    MWCLattice<Int, Dbl>* mwclat = 0;
    while (!timer.timeOver(timeLimit)) {
      if (type == MRG) {
        mrglat = nextGenerator(mrglat);
        Dbl merit(test(*mrglat));
        if (merit > bestMerit) {
          addLattice(mrglat, merit);
        }
      } else if (type == MWC) {
        mwclat = nextGenerator(mwclat);
        Dbl merit(test(*mwclat));
        if (merit > bestMerit) {
          addLattice(mwclat, merit);
        }
      } else if (type == MMRG) {
        mmrglat = nextGenerator(mmrglat);
        Dbl merit(test(*mmrglat));
        if (merit > bestMerit) {
          addLattice(mmrglat, merit);
        }
      }
      num_gen++;
    }
  }

  // Reads parameters from config file.
  bool readConfigFile(int argc, char **argv) {
    ParamReaderExt<Int, Dbl> reader(argv[1]);
    reader.getLines();
    int power;
    int ln = 0;
    reader.readGenType(type, ln++, 0);
    if (type == MRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
    } else if (type == MWC) {
      reader.readInt(power, ln++, 0);
      b = NTL::power2_ZZ(power);
      reader.readLong(exponent, ln++, 0);
    } else if (type == MMRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
    }
    reader.readInt(minDim, ln, 0);
    reader.readInt(maxDim, ln++, 1);
    reader.readDouble(timeLimit, ln++, 0);
    return true;
  }

}

int main (int argc, char **argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return -1;
  }
  readConfigFile(argc, argv);
  timer.init();
  testGenerators();
  printResults();
  delete bestLattice;
  //delete mod;
  return 0;
}

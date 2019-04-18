/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "SeekRe.h"

using namespace LatMRG;
using namespace LatticeTester::Random;

namespace {
  // Program global objects
  Chrono timer; // program timer
  LatticeTester::Normalizer<Dbl>* norma; // if a normalizer is used
  //IntLattice<Int, Int, Dbl, Dbl>* bestLattice = NULL; // the best lattice
  //Dbl bestMerit = Dbl(0); // best merit found
  long num_gen = 0;
  // For period tests
  DecompType decompm1 = DECOMP, decompr = DECOMP;
  std::string filem1, filer;
  TestList* bestLattice;

  // Data file parameters
  GenType type; // This one is experimental
  int minDim, maxDim;
  LatticeTester::NormaType normaType;
  double timeLimit;
  int max_gen;

  // MRG specific parameters
  MRGComponent<Int>* mrg;

  // MWC Specific parameters
  Int b; // modulo of MWC recurence

  // Shared components names
  std::int64_t order; // k for MRG and MMRG
  Int modulo; // m for MRG and MMRG
  // modulo is basis^exponent+rest
  std::int64_t basis;
  std::int64_t exponent;
  std::int64_t rest;
  bool period; // Period is full if this is true

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
    // Setting up two vectors. MRGComponent and MRGLattice do not use the same
    // vector format
    IntVec A;
    A.SetLength(order+1);
    NTL::clear(A);
    for (long i = 0; i<order; i++) A[i+1] = randInt(Int(0), modulo);
    // The program will not run the maxPeriod function if it is not wanted with
    // this condition
    while ((A[order] == 0) || (period && !mrg->maxPeriod(A))) {
      for (long i = 0; i<order; i++) A[i+1] = randInt(Int(0), modulo);
    }
    if (lattice) delete lattice;
    return new MRGLattice<Int, Dbl>(modulo, A, maxDim, order, FULL);
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

  void printResults() {
    std::cout << "Number of generators tested: " << num_gen << std::endl;
    auto list = bestLattice->getList();
    for (auto it = list.begin(); it!= list.end(); it++) {
      std::cout << "Message: " << (*it).getLattice() << "\nMerit: "
        << (*it).getMerit() << std::endl;
    }
    std::cout << "CPU time: " << timer.toString() << std::endl;
  }

  /**
   * Tests the generator via spectral test.
   * */
  DblVec test(IntLattice<Int, Int, Dbl, Dbl> & lattice) {
    norma = lattice.getNormalizer(normaType, 0, true);
    //Dbl merit = Dbl(1);
    DblVec results;
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
      //std::cout << "Pre-Normalization: " << NTL::sqrt(tmp) << std::endl;
      // Normalization
      //std::cout << "Bound: " << norma->getBound(i) << std::endl;
      tmp = NTL::sqrt(tmp)/norma->getBound(i);
      if (tmp > 1) tmp = Dbl(1)/tmp;
      //std::cout << "Value: " << tmp << std::endl;
      results.append(tmp);
      //merit = (tmp < merit) ? tmp : merit;
    }

    delete norma;
    return results;
  }

  // This is the main program loop. This loop searches for the next generator
  // and launches tests on it.
  void testGenerators() {
    if (type == MRG) {
      MRGLattice<Int, Dbl>* mrglat = 0;
      while (!timer.timeOver(timeLimit)) {
        mrglat = nextGenerator(mrglat);
        Test the_test("mrglattice", test(*mrglat));
        bestLattice->add(the_test);
        num_gen++;
      }
    } else if (type == MWC) {
      MWCLattice<Int, Dbl>* mwclat = 0;
      while (!timer.timeOver(timeLimit)) {
        mwclat = nextGenerator(mwclat);
        Test the_test("mwclattice", test(*mwclat));
        bestLattice->add(the_test);
        num_gen++;
      }
    } else if (type == MMRG) {
      MMRGLattice<Int, Dbl>* mmrglat = 0;
      while (!timer.timeOver(timeLimit)) {
        mmrglat = nextGenerator(mmrglat);
        Test the_test("mmrglattice", test(*mmrglat));
        bestLattice->add(the_test);
        num_gen++;
      }
    }
  }

  // Reads parameters from config file.
  bool readConfigFile(int argc, char **argv) {
    ParamReaderExt<Int, Dbl> reader(argv[1]);
    reader.getLines();
    int power;
    int ln = 0;
    reader.readGenType(type, ln++, 0);
    reader.readInt(minDim, ln, 0);
    reader.readInt(maxDim, ln++, 1);
    reader.readNormaType(normaType, ln++, 0);
    reader.readInt(max_gen, ln++, 0);
    reader.readDouble(timeLimit, ln++, 0);
    if (type == MRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
      minDim = minDim<order ? order : minDim;
      maxDim = maxDim>minDim ? maxDim : minDim;
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
      ln++;
    } else if (type == MWC) {
      reader.readInt(power, ln++, 0);
      b = NTL::power2_ZZ(power);
      reader.readLong(exponent, ln++, 0);
    } else if (type == MMRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
    }
    return true;
  }

}

int main (int argc, char **argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return -1;
  }
  // Initializing values
  srand(time(NULL));
  filem1 = "./tempm1" + std::to_string(rand());
  filer = "./tempr" + std::to_string(rand());
  readConfigFile(argc, argv);
  // Dynamically allocated objects
  mrg = new MRGComponent<Int>(modulo, order, decompm1, filem1.c_str(), decompr, filer.c_str());
  bestLattice = new TestList(max_gen);
  //else if (type == MWC) bestLattice = new TestList<MWCLattice>(max_gen);
  //else if (type == MMRG) bestLattice = new TestList<MMRGLattice>(max_gen);
  timer.init();
  // Launching the tests
  testGenerators();
  printResults();
  delete bestLattice;
  delete mrg;
  //delete mod;
  return 0;
}

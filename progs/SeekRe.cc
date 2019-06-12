/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#define LATMRG_SEEK
#include "Exec.h"

// Number of generators to generate before testing time in nextGenerator
#define DELAY 1000

using namespace LatMRG;
using namespace LatticeTester::Random;

namespace {
  // Program global objects
  Chrono timer; // program timer
  Normalizer<Dbl>* norma; // if a normalizer is used
  long num_gen = 0;
  // For period tests
  DecompType decompm1 = DECOMP, decompr = DECOMP;
  std::string filem1, filer;
  TestList* bestLattice;

  std::string delim = "\n========================================"
    "========================================\n\n";

  // Data file parameters
  GenType type;
  // Type of figure of merit
  LatticeTester::NormaType normaType = LatticeTester::NONE;
  LatticeTester::CriterionType criterion;
  LatticeTester::PreReductionType reduction;
  bool use_dual, best;
  // Projections
  int numProj;
  int minDim, maxDim;
  std::vector<std::size_t> projDim;
  Projections* proj;

  double timeLimit;
  int max_gen;

  // MRG specific parameters
  MRGComponent<Int>* mrg;
  int* coeff = NULL;

  // MWC Specific parameters
  Int b; // modulo of MWC recurence

  // Shared components names
  std::int64_t order; // k for MRG and MMRG
  Int modulo; // m for MRG and MMRG
  // modulo is basis^exponent+rest
  std::int64_t basis;
  std::int64_t exponent;
  std::int64_t rest;
  std::string construction;
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
    int delay = 0;
    // The program will not run the maxPeriod function if it is not wanted with
    // this condition
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<order; i++) A[i+1] = coeff[i] * randInt(Int(0), modulo);
      delay++;
    } while ((A[order] == 0) || (period && !mrg->maxPeriod(A)));
    if (lattice) delete lattice;
    return new MRGLattice<Int, Dbl>(modulo, A, maxDim, order, FULL);
  }

  MRGLattice<Int, Dbl>* nextGeneratorPow2(MRGLattice<Int, Dbl>* lattice) {
    // Setting up two vectors. MRGComponent and MRGLattice do not use the same
    // vector format
    IntVec A;
    A.SetLength(order+1);
    NTL::clear(A);
    int coefficients[2*order];
    int sign;
    int delay = 0;
    // The program will not run the maxPeriod function if it is not wanted with
    // this condition
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<order; i++) {
        if (coeff[2*i] < 0) {
          // This is a placeholder value for a zero coefficient
          coefficients[2*i] = coefficients[2*i+1] = 2004012;
          A[i+1] = 0;
          continue;
        }
        coefficients[2*i] = randInt(0, coeff[2*i]);
        sign = randInt(0,1);
        {
          Int tmp;
          NTL::power2(tmp, coefficients[2*i]);
          A[i+1] = Int(sign?1:-1) * tmp;
        }
        coefficients[2*i] ^= sign<<30;
        if (!(coeff[2*i+1] < 0)) {
          coefficients[2*i+1] = randInt(0, coeff[2*i+1]);
          sign = randInt(0,1);
          Int tmp;
          NTL::power2(tmp, coefficients[2*i+1]);
          A[i+1] += Int(sign?1:-1) * tmp;
          coefficients[2*i+1] ^= (sign<<30);
        }
        else coefficients[2*i+1] = 2004012;
      }
      delay++;
    } while ((A[order] == 0) || (period && !mrg->maxPeriod(A)));
    if (lattice) delete lattice;
    MRGLattice<Int, Dbl>* lat = new MRGLattice<Int, Dbl>(modulo, A, maxDim, order, FULL);
    lat->setPower2(coefficients);
    return lat;
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
    int delay = 0;
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<order-1; i++) {
        A[i][i+1] = Int(1);
      }
      for (long i = 0; i<order; i++) {
        A[order-1][i] = randInt(Int(0), modulo);
      }
      delay++;
    } while ((NTL::determinant(A) == 0) || (period && !mrg->maxPeriod(A)));
    // Correcting the matrix to a full matrix
    for (int i = 0; i<order; i++) A *= A;
    for (int i = 0; i<order; i++)
      for (int j = 0; j<order; j++) 
        A[i][j] = A[i][j]%modulo;
    if (lattice) delete lattice;
    return new MMRGLattice<Int, Dbl>(modulo, A, maxDim, order);
  }

  /*
   * These next function add the tested lattices to the list of the best ones.
   * This only add the lattices that are good enough.
   * */
  void printResults() {
    std::cout << "\nSeekRe: A search program for Random Number Generators\n";
    std::cout << delim;
    std::cout << "Bellow are the results of a search for random number generators:\n";
    std::cout << "Generator type: " << toStringGen(type) << "\n";
    if (type == MRG) {
      std::cout << "With modulus:   m = " << modulo << " = " << basis << "^"
        << exponent << (rest>0?"+":"") << (rest!=0?std::to_string(rest):"") << "\n";
      std::cout << "Of order:       k = " << order << "\n";
    } else if (type == MWC) {
    } else if (type == MMRG) {
    }
    std::cout << "And " << (period?"full":"any") << " period length\n";
    std::cout << "The test was:\n" << (best?"Best":"Worst") << " generators "
      "ranked by ";
    if(criterion == LatticeTester::SPECTRAL) std::cout << "minimal " 
      << (normaType==LatticeTester::NONE?"inverse":"normalized")
      << " shortest non-zero vector length (Spectral test)\n";
    else if (criterion == LatticeTester::LENGTH) std::cout << "minimal"
      << " shortest non-zero vector length\n";
    else if (criterion == LatticeTester::BEYER) std::cout << "their Beyer quotient\n";
    if (normaType != LatticeTester::NONE) {
      std::cout << "Normalizer used: "
        << LatticeTester::toStringNorma(normaType) << "\n";
    }
    std::cout << "On dimensions and projections:\n";
    std::cout << proj->toString();
    std::cout << delim;
    std::cout << "Allowed running time: " << timeLimit << "s.\n";
    std::cout << "Actual CPU time: " << timer.toString() << "\n";
    std::cout << "Number of generators kept: " << max_gen << "\n";
    std::cout << "Number of generators tested: " << num_gen << "\n\n";
    std::cout << "Retained generators (from best to worst):\n";
    auto list = bestLattice->getList();
    for (auto it = list.begin(); it!= list.end(); it++) {
      std::cout << delim;
      if (type == MRG) {
        std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
      } else if (type == MWC) {}
      else if (type == MMRG) {
        std::cout << "Matrix:\n" << (*it).getLattice() << "\n";
      }
      std::cout << "Merit: " << (*it).getMerit() << "\n";
    }
  }

  /**
   * Tests the generator via spectral test.
   * */
  DblVec test(IntLattice<Int, Int, Dbl, Dbl> & lattice) {
    norma = lattice.getNormalizer(normaType, 0, true);
    //Dbl merit = Dbl(1);
    DblVec results(0);
    lattice.buildBasis(minDim);
    for (int i = minDim; i <= maxDim; i++) {
      // Changing to the dual
      if (use_dual) lattice.dualize();
      // Reducing the lattice
      if (reduction == LatticeTester::FULL)
        reduceFull(lattice);
      else if (reduction == LatticeTester::LLL)
        reduceLLL(lattice);
      else if (reduction == LatticeTester::BKZ)
        reduceBKZ(lattice);
      else if (reduction == LatticeTester::NOPRERED)
        reduceMink(lattice);
      // Computing the merit of the lattice
      Dbl tmp;
      if (criterion == LatticeTester::LENGTH) tmp = meritL(lattice, norma);
      if (criterion == LatticeTester::SPECTRAL) tmp = meritS(lattice, norma);
      if (criterion == LatticeTester::BEYER) tmp = meritB(lattice, norma);
      results.append(tmp);
      // Rejecting lattices that won't make it
      if (tmp < bestLattice->getMerit()) {
        results[0] = 1-best;
        return results;
      }
      // Changing back to the primal and increasing the dimension
      if (use_dual) lattice.dualize();
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
        if (use_dual) proj_lat.dualize();
        // Reduction
        if (reduction == LatticeTester::FULL)
          reduceFull(proj_lat);
        else if (reduction == LatticeTester::LLL)
          reduceLLL(proj_lat);
        else if (reduction == LatticeTester::BKZ)
          reduceBKZ(proj_lat);
        else if (reduction == LatticeTester::NOPRERED)
          reduceMink(proj_lat);
        // Figure of merit
        Dbl tmp;
        if (criterion == LatticeTester::LENGTH) tmp = meritL(proj_lat, norma);
        if (criterion == LatticeTester::SPECTRAL) tmp = meritS(proj_lat, norma);
        if (criterion == LatticeTester::BEYER) tmp = meritB(proj_lat, norma);
        results.append(tmp);
        // Rejecting lattices that won't make it
        if (tmp < bestLattice->getMerit()) {
          results[0] = 1-best;
          return results;
        }
      }
    }

    delete norma;
    return results;
  }

  int print_progress(int old) {
    int per_80 = timer.val(Chrono::SEC)/timeLimit * 80;
    if (per_80 > 80) per_80 = 80;
    if (per_80 < 0) per_80 = 0;
    // We do not print for no reason as this slows the program a lot.
    if (per_80 <= old) return old;
    std::cout << "[";
    for (int i = 0; i < per_80; i++) std::cout << "#";
    for (int i = per_80; i < 80; i++) std::cout << " ";
    std::cout << "] ";
    std::cout << std::setw(2) << int(per_80/80.0*100) << " %\r" << std::flush;
    return per_80;
  }

  // This is the main program loop. This loop searches for the next generator
  // and launches tests on it.
  void testGenerators() {
    std::cout << "Program progress:\n";
    int old = print_progress(-1);
    if (type == MRG) {
      MRGLattice<Int, Dbl>* mrglat = 0;
      while (!timer.timeOver(timeLimit)) {
        //std::cout << "1"<<std::endl;
        if (construction == "POW2") mrglat = nextGeneratorPow2(mrglat);
        else if (construction == "RANDOM") mrglat = nextGenerator(mrglat);
        //std::cout << "2" << std::endl;
        if (construction == "POW2") mrglat = nextGeneratorPow2(mrglat);
        if (mrglat == NULL) continue;
        Test the_test(mrglat->toString(), test(*mrglat));
        bestLattice->add(the_test);
        num_gen++;
        old = print_progress(old);
      }
    } else if (type == MWC) {
      MWCLattice<Int, Dbl>* mwclat = 0;
      while (!timer.timeOver(timeLimit)) {
        mwclat = nextGenerator(mwclat);
        if (mwclat == NULL) continue;
        Test the_test("mwclattice", test(*mwclat));
        bestLattice->add(the_test);
        num_gen++;
        old = print_progress(old);
      }
    } else if (type == MMRG) {
      MMRGLattice<Int, Dbl>* mmrglat = 0;
      while (!timer.timeOver(timeLimit)) {
        mmrglat = nextGenerator(mmrglat);
        if (mmrglat == NULL) continue;
        Test the_test(mmrglat->toStringGeneratorMatrix(), test(*mmrglat));
        bestLattice->add(the_test);
        num_gen++;
        old = print_progress(old);
      }
    }
  }

  // Reads parameters from config file.
  bool readConfigFile(int argc, char **argv) {
    ParamReaderExt<Int, Dbl> reader(argv[1]);
    reader.getLines();
    int power;
    int ln = 0;
    // Reading global problem parameters
    reader.readGenType(type, ln++, 0);
    reader.readCriterionType(criterion, ln, 0);
    if (criterion == LatticeTester::SPECTRAL) {
      reader.readBool(use_dual, ln, 1);
      reader.readPreRed(reduction, ln, 2);
      reader.readNormaType(normaType, ln++, 3);
    } else if (criterion == LatticeTester::LENGTH) {
      reader.readBool(use_dual, ln, 1);
      reader.readPreRed(reduction, ln++, 2);
    } else if (criterion == LatticeTester::BEYER) {
      reduction = LatticeTester::NOPRERED;
      ln++;
    }
    reader.readBool(best, ln++, 0);
    // This code corrects the projections so that it always builds a valid
    // projection specification
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
    if (numProj > 1) ln++;
    reader.readInt(max_gen, ln++, 0);
    reader.readDouble(timeLimit, ln++, 0);
    if (type == MRG) {
      reader.readNumber3(modulo, basis, exponent, rest, ln++, 0);
      reader.readLong(order, ln++, 0);
      // Reading the construction method
      reader.readString(construction, ln, 0);
      if (construction == "RANDOM") {
        coeff = new int[order];
        for (unsigned int i = 1; i < order; i++) reader.readInt(coeff[i-1], ln, i);
        coeff[order-1] = 1;
        ln++;
      } else if (construction == "POW2") {
        coeff = new int[2 * order];
        for (unsigned int i = 1; i < 2*order-1; i++)
          reader.readInt(coeff[i-1], ln, i);
        coeff[2*order-1] = exponent-1;
        coeff[2*order-2] = exponent-1;
        ln++;
      }
      reader.readBool(period, ln, 0);
      // Making sure that minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      minDim = (!period||minDim>order) ? minDim : (order+1);
      maxDim = maxDim>minDim ? maxDim : minDim;
      for (unsigned int i = 0; i < projDim.size(); i++) {
        projDim[i] = (projDim[i]<(unsigned)(minDim-1))?(unsigned)(minDim-1):projDim[i];
      }
      if (period) {
        // Using default parameters for period or not
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
      reader.readBool(period, ln, 0);
      // Making sure that minDim is big enough to provide usefull tests (if
      // full period is required) this changes other dimensions accordingly
      minDim = (!period||minDim>order) ? minDim : (order+1);
      maxDim = maxDim>minDim ? maxDim : minDim;
      for (unsigned int i = 0; i < projDim.size(); i++) {
        projDim[i] = (projDim[i]<(unsigned)(minDim-1))?(unsigned)(minDim-1):projDim[i];
      }
      if (period) {
        // Using default parameters for period or not
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
  bestLattice = new TestList(max_gen, best);
  proj = new Projections(numProj, minDim, projDim);
  //else if (type == MWC) bestLattice = new TestList<MWCLattice>(max_gen);
  //else if (type == MMRG) bestLattice = new TestList<MMRGLattice>(max_gen);
  timer.init();
  // Launching the tests
  testGenerators();
  printResults();
  delete proj;
  delete bestLattice;
  delete mrg;
  if (coeff) delete[] coeff;
  return 0;
}

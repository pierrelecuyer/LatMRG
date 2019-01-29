/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "latticetester/NormaBestLat.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/ParamReaderExt.h"

using namespace LatMRG;

namespace {
  typedef NTL::ZZ Int;
  typedef NTL::RR Dbl;

  Chrono timer;
  LatticeTester::Normalizer<Dbl>* norma;
  MWCLattice<Int, Dbl>* bestLattice = NULL;
  Dbl bestMerit = Dbl(0);
  long num_gen = 0;

  Int b = NTL::power2_ZZ(64);
  int exponent;
  int minDim, maxDim;
  double timeLimit;

  /**
   * A small class to search for modulus for MWC generators.
   * */
  class Modulus {
    public:

      /**
       * Search will begin at `2^e` and will increase of decrease depending on
       * `increase`.
       * */
      Modulus(long e, bool increase) {
        m = (Int(1)<<e) - 1;
        this->increase = increase;
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
  } *mod;

  /**
   * Tests the generator via spectral test.
   * */
  Dbl test(MWCLattice<Int, Dbl> & lattice) {
    norma = lattice.getNormalizer(LatticeTester::BESTLAT, 0, true);
    Dbl merit = Dbl(1);
    for (int i = minDim; i <= maxDim; i++){
      std::cout << "i: " << i << std::endl;
      // Building the basis
      lattice.buildBasis(i);
      lattice.dualize();
      // Reducing the lattice
      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lattice);
      red.redBKZ(0.999999, 10, LatticeTester::EXACT, lattice.getDim());
      red.shortestVector(lattice.getNorm());
      // Computing shortest vector length and spectral test
      NTL::vector<Int> shortest(lattice.getBasis()[0]);
      Dbl tmp;
      LatticeTester::ProdScal<Int>(shortest, shortest, i, tmp);
      std::cout << "Pre-Normalization: " << NTL::sqrt(tmp) << std::endl;
      // Normalization
      std::cout << "Bound: " << norma->getBound(i) << std::endl;
      tmp = NTL::sqrt(tmp)/norma->getBound(i);
      std::cout << "Normalized value: " << tmp << std::endl;
      merit = (tmp < merit) ? tmp : merit;
    }

    delete norma;
    return merit;
  }


  // This just instanciates number MWC generators with order k and mod b.
  void testGenerators() {
    while (!timer.timeOver(timeLimit)) {
      MWCLattice<Int, Dbl> lattice(b, mod->next());
      Dbl merit(test(lattice));
      if (merit > bestMerit) {
        bestMerit = merit;
        if(bestLattice) delete bestLattice;
        bestLattice = new MWCLattice<Int, Dbl>(lattice);
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
    reader.readInt(power, ln++, 0);
    b = NTL::power2_ZZ(power);
    reader.readInt(exponent, ln++, 0);
    mod = new Modulus(exponent, true);
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
  std::cout << "Number of generators tested: " << num_gen << std::endl;
  std::cout << "Best merit: " << bestMerit << std::endl;
  std::cout << "Best lattice LCG coefficient: " << bestLattice->getCoef() << std::endl;
  std::cout << "Best lattice modulus: " << bestLattice->getModulo() << std::endl;
  std::cout << "CPU time: " << timer.toString() << std::endl;
  delete bestLattice;
  delete mod;
  return 0;
}

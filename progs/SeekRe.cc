/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "latmrg/MWCLattice.h"

using namespace LatMRG;

namespace {
  NTL::ZZ b = NTL::power2_ZZ(64);
  long exponent = 250;

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
        m = NTL::ZZ(1)<<e;
        this->increase = increase;
      }

      /**
       * This function finds the next value for `m`. This can return 1 if
       * `increase` is `false` and `m` gets to 1.
       * */
      NTL::ZZ next() {
        while (m > 1) {
          nextM();
          LatticeTester::PrimeType status = LatticeTester::IntFactor<NTL::ZZ>::isPrime (m, KTRIALS);
          if (status == LatticeTester::PRIME || status == LatticeTester::PROB_PRIME) {
            NTL::ZZ m1 = (m - 1)/2;
            status = LatticeTester::IntFactor<NTL::ZZ>::isPrime (m1, KTRIALS);
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
        return NTL::ZZ(-1);
      }

    private:

      static const long KTRIALS = 200;

      /**
       * Last valid modulus found (or 2^e).
       * */
      NTL::ZZ m;

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
  } mod(exponent, true);
  /**
   * Reads the configuration file for this program.
   * */
  bool readConfigFile(int argc, char **argv) {
    return true;
  }

  // This just instanciates number MWC generators with order k and mod b.
  void spawnGenerators(int number) {
    for (long i = 0; i<number; i++) {
      MWCLattice<NTL::ZZ, double> lattice(b, mod.next());
      std::cout << "Modulus: " << lattice.getCoef() << std::endl;
    }
  }
}

int main (int argc, char **argv)
{
  readConfigFile(argc, argv);
  spawnGenerators(10);
  //if (readConfigFile (argc, argv))
  //  return -1;
  //config.write ();
  // Init ();
  // timer.init ();
  // Zone<MScal>::initFrontieres (config);
  // InitZones ();
  // if (config.readGenFile) {
  //   TestGen();
  // } else {
  //   SeekGen (0);
  // }
  // SortBestGen ();
  // if (config.outputType == RES) {
  //   PrintResults ();
  // } else if (config.outputType == GEN) {
  //   PrintGen();
  // } else {
  //   std::cout << "\nCould not print results as requested, tried to output using
  //     standard method.\n";
  //   PrintResults ();
  // }
  // Finalize ();
  return 0;
}

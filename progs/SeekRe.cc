/**
 * This program rewrites parts of SeekMain to be easier to understand and use
 * better design. For now this is the file in which I implement MWC gen searches.
 * */
#include "latmrg/MWCLattice.h"

using namespace LatMRG;

namespace {
  int k = 5;
  NTL::ZZ b = NTL::power2_ZZ(64);
  /**
   * Reads the configuration file for this program.
   * */
  bool readConfigFile(int argc, char **argv) {
    return true;
  }

  // This just instanciates number MWC generators with order k and mod b.
  void spawnGenerators(int number) {
    NTL::vector<NTL::ZZ> coeff;
    coeff.SetLength(k+1);
    coeff[0] = 3;
    coeff[1] = 4;
    coeff[2] = 7;
    coeff[3] = 11;
    coeff[4] = 13;
    coeff[5] = 17;
    for (long i = 0; i<number; i++) {
      for(int j = 0; j <= k; j++) {
        coeff[j] = NTL::ZZ(LatticeTester::RandBits(63));
        std::cout << coeff[j] << std::endl;
      }
      MWCLattice<NTL::ZZ, double> lattice(b, coeff, NTL::ZZ(1), k);
      std::cout << lattice.toStringCoef() << std::endl;
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

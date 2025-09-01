/*
 * TestMRG32k3a.cc
 */

// This defines the Int type. We must recompile to change it.
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Chrono.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latmrg/MRGLattice.h"
// #include "latmrg/MRGLatticeLac.h"

using namespace LatticeTester;
using namespace LatMRG;

int main() {
   // Note how large ZZ integers can be initialized in NTL.
   Int m = to_ZZ("18446645023178547541");
   IntVec aa;
   aa.SetLength(4);
   aa[1] = to_ZZ("18169668471252892557");
   aa[2] = to_ZZ("3186860506199273833");
   aa[3] = to_ZZ("8738613264398222622");

   int64_t maxdim(50);  // Maximum dimension of the lattice.  ***

   NTL::Vec <int64_t> t; // The t-vector for the FOM.
   t.SetLength(3);
   t[0] = maxdim;  t[1] = 0;  t[2] = 0;
   MRGLattice<Int, Real> mrgc = MRGLattice<Int, Real>(m, aa, maxdim, 0, maxdim);
   WeightsUniform weights(1.0);
   NormaRogers normaDual(-log(m), 3, maxdim);  // Normalization for m-dual: Roger's bounds as in rLEC99b.
   ReducerBB<Int, Real> red(maxdim);           // Single ReducerBB.
   FigureOfMeritDualM<Int, Real> fom(t, weights, normaDual, &red, true); // FoM for dual lattice.
   fom.setVerbosity(3);
   std::cout << "\nResults from TestMRG32k3a.cc \n";

   // We first compute the FOM in successive dimensions up to maxdim, with the default BKZ + BB.
   // It uses delta = 0.99999 and blocksize = 10.
   fom.setBKZ(0.99999, 10);
   std::cout << "\n=============================================\n";
   std::cout << "MRG32k3a, with BKZ+BB with default parameters \n";
   fom.computeMeritSucc(mrgc);

   // Again BKZ + BB, but with delta = 0.99999 and blocksize k = 20.
   fom.setBKZ(0.99999, 20);
   std::cout << "\n=============================================\n";
   std::cout << "MRG32k3a, with BKZ+BB with delta = 0.99999 and block size = 20 \n";
   fom.computeMeritSucc(mrgc);

   // Again BKZ + BB, but with delta = 0.99999999 and blocksize k = 50.
   fom.setBKZ(0.99999999, 50);
   std::cout << "\n=============================================\n";
   std::cout << "MRG32k3a, with BKZ+BB with delta = 0.99999999 and block size = 50 \n";
   fom.computeMeritSucc(mrgc);

   // Here we apply only LLL with delta = 0.9, no BB.
   fom.setBB(false);
   fom.setBKZ(0.0);
   fom.setLLL(0.9);
   std::cout << "\n===========================================\n";
   std::cout << "MRG32k3a, only LLL with delta = 0.9, no BB \n";
   fom.computeMeritSucc(mrgc);

   // Here we apply only the default BKZ, no BB.
   fom.setBKZ(0.99999);
   std::cout << "\n===========================================\n";
   std::cout << "MRG32k3a, only BKZ with block size k = 10, no BB \n";
   fom.computeMeritSucc(mrgc);

   // Here we apply only BKZ, but with delta = 0.99999999 and blocksize k = 50.
   fom.setBKZ(0.99999999, 50);
   std::cout << "\n===========================================\n";
   std::cout << "MRG32k3a, only BKZ with delta = 0.99999999 and k = 50, no BB \n";
   fom.computeMeritSucc(mrgc);

   // Here we rebuild the basis from scratch each time we increase the dimension,
   // with the default BKZ + BB.
   fom.setLLL(0.0);
   fom.setBKZ(0.99999, 10);
   fom.setBB(true);
   std::cout << "\n===========================================\n";
   std::cout << "MRG32k3a, default BKZ+BB, rebuild basis at each step \n";
   fom.computeMeritSuccRebuild(mrgc);

   // Rebuild basis from scratch and use only LLL with delta = 0.9.
   fom.setLLL(0.9);
   fom.setBKZ(0.0);
   fom.setBB(false);
   std::cout << "\n===========================================\n";
   std::cout << "MRG32k3a, only LLL with delta = 0.9, rebuild basis at each step \n";
   fom.computeMeritSuccRebuild(mrgc);

   return 0;
}



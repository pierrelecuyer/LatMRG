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
#include "latmrg/MRGLatticeLac.h"

using namespace LatticeTester;
using namespace LatMRG;

int main() {
   Int m = to_ZZ("18446645023178547541");
   IntVec aa;
   aa.SetLength(4);
   aa[1] = to_ZZ("18169668471252892557");
   aa[2] = to_ZZ("3186860506199273833");
   aa[3] = to_ZZ("8738613264398222622");
   Chrono timer;

   int64_t maxdim(45);  // Maximum dimension of the lattice
   NTL::Vec <int64_t> t; // The t-vector for the FOM.
   t.SetLength(3);
   t[0] = maxdim;  t[1] = 0;  t[2] = 0;
   double merit, minmerit;
   MRGLattice<Int, Real> mrgc = MRGLattice<Int, Real>(m, aa, maxdim);
   WeightsUniform weights(1.0);
   NormaRogers normaDual(-log(m), 3, maxdim);  // Normalization for m-dual: Roger's bounds as in rLEC99b.
   ReducerBB<Int, Real> red(maxdim);   // Single ReducerBB with internal lattice `lat`.
   FigureOfMeritDualM<Int, Real> fom(t, weights, normaDual, &red, true); // FoM for dual lattice.
   fom.setVerbosity(3);
   std::cout << "\nResults from TestMRG32k3a.cc \n\n";

   // We first compute the FOM in successive dimensions up to maxdim, with the default BKZ + BB.
   // It uses delta = 0.99999 and blocksize = 10.
   std::cout << "=============================================\n";
   std::cout << "MRG32k3a, with BKZ+BB with default parameters \n";
   timer.init();
   mrgc.buildDualBasis(4);
   minmerit = fom.computeMeritSucc(mrgc);
   std::cout << "minmerit = " << minmerit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Here we apply only the default BKZ, no BB.
   fom.setBB(false);
   std::cout << "===========================================\n";
   std::cout << "MRG32k3a, only BKZ, no BB \n";
   timer.init();
   mrgc.buildDualBasis(4);
   minmerit = fom.computeMeritSucc(mrgc);
   std::cout << "minmerit = " << minmerit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Then only LLL with delta = 0.9.
   fom.setBKZ(0.0);
   fom.setLLL(0.9);
   std::cout << "===========================================\n";
   std::cout << "MRG32k3a, only LLL with delta = 0.9, no BB \n";
   timer.init();
   mrgc.buildDualBasis(4);
   minmerit = fom.computeMeritSucc(mrgc);
   std::cout << "minmerit = " << minmerit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Here we rebuild the basis from scratch each time we increase the dimension,
   // with the default BKZ + BB.
   int64_t j;
   fom.setLLL(0.99999);
   fom.setBKZ(0.99999, 10);
   fom.setBB(true);
   Coordinates coord;
   minmerit = DBL_MAX;
   std::cout << "===========================================\n";
   std::cout << "MRG32k3a, BKZ+BB, rebuild basis at each step \n";
   std::cout << "coordinates      sqlen       merit           minmerit \n";
   timer.init();
   for (j = 1; j <= 3; j++) coord.insert(j);
   for (j = 4; j <= maxdim; j++) {
      coord.insert(j);
      mrgc.buildDualBasis(j);
      mrgc.dualize();
      merit = fom.computeMeritOneProj(mrgc, coord, minmerit);
      minmerit = min(merit, minmerit);
   }
   std::cout << "minmerit = " << minmerit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Rebuild basis from scratch and use only LLL with delta = 0.9.
   fom.setLLL(0.9);
   fom.setBKZ(0.0);
   fom.setBB(false);
   minmerit = DBL_MAX;
   coord.clear();
   std::cout << "===========================================\n";
   std::cout << "MRG32k3a, only LLL with delta = 0.9, rebuild basis at each step \n";
   std::cout << "coordinates      sqlen       merit           minmerit \n";
   timer.init();
   for (j = 1; j <= 3; j++) coord.insert(j);
   for (j = 4; j <= maxdim; j++) {
      coord.insert(j);
      mrgc.buildDualBasis(j);
      mrgc.dualize();
      merit = fom.computeMeritOneProj(mrgc, coord, minmerit);
      minmerit = min(merit, minmerit);
   }
   std::cout << "minmerit = " << minmerit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   return 0;
}



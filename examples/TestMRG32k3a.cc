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

   int64_t maxdim(46);  // Maximum dimension of the lattice
   NTL::Vec <int64_t> t; // The t-vector for the FOM.
   t.SetLength(3);
   t[0] = maxdim;  t[1] = 0;  t[2] = 0;
   MRGLattice<Int, Real> mrgc = MRGLattice<Int, Real>(m, aa, maxdim);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), 3, maxdim);  // Factors computed for dual.
   ReducerBB<Int, Real> red(maxdim);   // Single ReducerBB with internal lattice `lat`.
   FigureOfMeritDualM<Int, Real> fom(t, weights, normaDual, &red, true); // FoM for dual lattice.
   fom.setVerbosity(3);

   // We compute the FOM in successive dimensions up to maxdim, without applying the BB.
   fom.setBB(false);
   timer.init();
   mrgc.buildDualBasis(4);
   double merit = fom.computeMeritSucc(mrgc);
   std::cout << "MRG32k3a, no BB, merit = " << merit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Here we also apply the BB, also in up to maxdim dimensions.
   fom.setBB(true);
   timer.init();
   mrgc.buildDualBasis(4);
   merit = fom.computeMeritSucc(mrgc);
   std::cout << "MRG32k3a, with BB, merit = " << merit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   // Here we compute the shortest dual vector only in maxdim dimensions.
   timer.init();
   int64_t newdim = maxdim;
   mrgc.buildDualBasis(newdim);
   mrgc.dualize();
   Coordinates coord;   // We have to do this due to some bad design...
   for (int64_t j = 1; j <= newdim; j++)
      coord.insert(j);
   merit = fom.computeMeritOneProj(mrgc, coord);
   std::cout << "MRG32k3a, direct with BB, merit = " << merit << "\n";
   std::cout << "CPU time: " << timer.toString() << "\n\n";

   return 0;
}



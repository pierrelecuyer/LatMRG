
/**
 * This example shows how projections sometimes have a smaller density than the full lattice.
 * We do this with a small MRG example or order 3 with `m = 13`, then with a larger MRG example
 * with `m` near @f$2^{63}@f$. For each case, we compute the FOM M_t for a given set of projections,
 * and we find that the projection over coordinates {1,3,4} has a smaller density because some
 * points project only each other.
 */
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>


#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"
#include "latmrg/MRGLattice.h"

using namespace LatticeTester;
// using namespace LatMRG;

// Tests an MRG of modulus m and multipliers a, for projections specified by t, with
// coordinate 1 always included.  The order k is the size of the vector a.
// We first examine the projections over successive coordinates, then the other ones.
// We do this for the primal, then for the m-dual.
// Finally, we examine more closely the projection over {1, 3, 4}.
template<typename Int, typename Real>
void testProjectionsMRG (const Int m, const NTL::Vec<Int> a, const NTL::Vec<int64_t> t) {
   
   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   int64_t order = a.length()-1;
   double merit;
   LatMRG::MRGLattice<Int, Real> lat(m, a, maxdim);
   WeightsUniform weights(1.0);
   ReducerBB<Int, Real> red(maxdim);   // Reducer created for up to maxdim dimensions.

   // Building full primal lattice and looking at projections over successive coordinates.
   std::cout << "===================================================\n";
   std::cout << "We build the lattice and look at projections over successive coordinates.\n";
   std::cout << "\nFigure of merit primal succ, with BB.\n";
   // We consider only the projections that contain coordinate 1.
   NormaBestLat norma(log(m), order, maxdim);  // Factors will be computed for primal.
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red, true);
   fom.setVerbosity(2);
   merit = fom.computeMeritSucc(lat);
   std::cout << "FOM value: " << merit << "\n\n";

   // Now looking at other primal projections, placed in lattice proj.
   std::cout << "===================================================\n";
   std::cout << "Then we look at other primal projections, over pairs and triples.\n";
   lat.buildBasis(maxdim);
   std::cout << "\nFigure of merit primal non-succ, with BB.\n";
   IntLattice<Int, Real> proj(m, 3);
   merit = fom.computeMeritNonSucc(lat, proj);
   std::cout << "FOM value: " << merit << "\n\n";

   // We now examine the m-dual lattice, first for successive coordinates.
   std::cout << "===================================================\n";
   std::cout << "We now examine the m-dual lattice, first for successive coordinates.\n";
   NormaBestLat normadual(-log(m), order, maxdim);  // Factors will be computed for m-dual.
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normadual, &red, true);
   fomdual.setVerbosity(2);
   std::cout << "\nFigure of merit for dual.\n";
   merit = fomdual.computeMeritSucc(lat);
   std::cout << "FOM value: " << merit << "\n\n";

   // Then the m-duals of the projections.
   std::cout << "===================================================\n";
   std::cout << "Then we look at the m-duals of the other projections.\n";
   std::cout << "\nFigure of merit dual non-succ, with BB.\n";
   merit = fomdual.computeMeritNonSucc(lat, proj);
   std::cout << "FOM value: " << merit << "\n\n";

   // A closer look at projection {1,3,4}
   std::cout << "===================================================\n";
   std::cout << "A closer look at the projection over coordinates {1,3,4}: \n";
   Coordinates coord({1, 3, 4});
   std::cout << "Full basis B before taking projection {1,3,4}: \n" << lat.getBasis() << "\n";
   lat.buildProjection(proj, coord);
   std::cout << "Basis for projection {1,3,4}: \n" << proj.getBasis() << "\n";
   fom.setVerbosity(0);
   fom.computeMeritOneProj(proj, coord);

   lat.buildProjectionDual(proj, coord);
   proj.dualize();
   std::cout << "Basis for m-dual projection {1,3,4}: \n" << proj.getBasis() << "\n";
   fomdual.setVerbosity(0);
   fomdual.computeMeritOneProj(proj, coord);
   
}


int main() {


   
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::Vec<int64_t> t; // The t-vector for the FOM.
   t.SetLength(3);
   t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 5;    // Then pairs and triples up to coord. 5.
   t[2] = 5;

   // This is the small example at the end of Section 9 of the guide.
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a small MRG example with m=13, k=3, a=(7,0,4). \n";
   NTL::ZZ m(13);
   NTL::Vec<NTL::ZZ> a; // Vector a has size 4, a[j] contains a_j.
   a.SetLength(4);
   a[1] = 7;
   a[2] = 0;
   a[3] = 4;
   testProjectionsMRG<NTL::ZZ, double> (m, a, t);

    // This is the MRG retained in Table VII of the LatMRG paper of L'Ecuyer and Couture (1997).
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   a[1] = 1145902849652723;
   a[2] = 0;
   a[3] = -1184153554609676;
   testProjectionsMRG<NTL::ZZ, double> (m, a, t);
   

   return 0;
}

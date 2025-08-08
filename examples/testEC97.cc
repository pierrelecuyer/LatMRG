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
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/ReducerStatic.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"

using namespace LatticeTester;

template<typename Int, typename Real>
static void FOMLattice(Int m, int64_t maxdim, const NTL::Vec<Int> a, const NTL::Vec<int64_t> t) {
  LatMRG::MRGLattice<NTL::ZZ, double> lat(m, a, maxdim);
  lat.setUsePolynomialBasis(true);
  ReducerBB<NTL::ZZ, double> m_red(lat); 
  WeightsUniform weights(1.0);
  NormaBestLat normaDual(-log(m), a.length()-1, maxdim, L2NORM);
  FigureOfMeritDualM<NTL::ZZ, double> fomdual(t, weights, normaDual, &m_red);
  fomdual.setVerbosity(2);
   IntLattice<NTL::ZZ, double> proj(m, t.length());
  fomdual.computeMerit(lat, proj);
}


int main() {

   int64_t maxdim = 8;   
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::Vec<int64_t> t; // The t-vector for the FOM.
   t.SetLength(2);
   t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 0;    // Then pairs and triples up to coord. 5.   
   NTL::ZZ m;
   NTL::Vec<NTL::ZZ> a; // Vector a has size 4, a[j] contains a_j.

   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 2147483647 and a = 45991\n";
   m = 2147483647;
   a.SetLength(2);
   a[1] = 45991;
   FOMLattice<NTL::ZZ, double>(m, maxdim, a, t);   
   std::cout << "\n The results are the same as in [EC97], Table 6!\n";


   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   a.SetLength(4);
   a[1] = 1145902849652723;
   a[2] = 0;
   a[3] = -1184153554609676;
  FOMLattice<NTL::ZZ, double>(m, maxdim, a, t);
  std::cout << "\n The results are the same as in [EC97], Table 7!\n";

   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 4294967296, a = 1099087573. \n";
   maxdim = 8;
   m = 4294967296;
   m = m / 4;
   a.SetLength(2);
   a[1] = 1099087573;  
   FOMLattice<NTL::ZZ, double>(m, maxdim, a, t);    
  std::cout << "\n The results are the same as in [EC97], Table 1!\n";

}

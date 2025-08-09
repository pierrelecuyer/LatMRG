/**
 * In this example, we show that the results of the new LatMRG version
 * are consistent with the examples reported in the the article [rLEC97c].
 * It is important that there are either no difference or, if they appear,
 * it needs to be shown that the calculations in the current version
 * are correct.
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

/**
 * This function calculates the FoM of the dual lattice of an MRG with modulus 'm' 
 * and vector of multipliers 'a'. The maximal dimension of the lattice is defined by 'maxdim' 
 * and the actual FoM is determined by the vector 't'. 
 */
template<typename Int, typename Real>
static void FoMMRGLattice(Int m, int64_t maxdim, const NTL::Vec<Int> a, const NTL::Vec<int64_t> t) {
   LatMRG::MRGLattice<Int, Real> lat(m, a, maxdim);
   lat.setUsePolynomialBasis(true);
   ReducerBB<Int, Real> m_red(lat); 
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), a.length()-1, maxdim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &m_red);
   fomdual.setVerbosity(2);
    IntLattice<Int, Real> proj(m, t.length());
   fomdual.computeMerit(lat, proj);
}

/**
 * This function does the same as FoMMRGLattice except that it deals with MRG lattices with
 * lacunary indices defined by 'lac'.
 */
template<typename Int, typename Real>
static void FoMMRGLatticeLac(Int m, int64_t maxdim, const NTL::Vec<Int> a, const NTL::Vec<int64_t> t, 
  NTL::Vec<Int> lac) {
   LatMRG::MRGLatticeLac<Int, Real> lat(m, a, maxdim, lac);
   ReducerBB<Int, Real> m_red(lat); 
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), a.length()-1, maxdim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &m_red);
   fomdual.setVerbosity(2);
   IntLattice<Int, Real> proj(m, t.length());
   fomdual.computeMerit(lat, proj);
}


int main() {

   int64_t maxdim = 8;   
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::Vec<int64_t> t; // The t-vector for the FoM.
   t.SetLength(2);
   t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 0;      
   NTL::ZZ m;
   NTL::Vec<NTL::ZZ> a; // Vector a has size d+1, if it contains d elements.
   
   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 2147483647 and a = 45991\n";
   m = 2147483647;
   a.SetLength(2);
   a[1] = 45991;
   FoMMRGLattice<NTL::ZZ, double>(m, maxdim, a, t);   
   std::cout << "\n The results are the same as in [rLEC97c], Table 6!\n";


   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   a.SetLength(4);
   a[1] = 1145902849652723;
   a[2] = 0;
   a[3] = -1184153554609676;
  FoMMRGLattice<NTL::ZZ, double>(m, maxdim, a, t);
  std::cout << "\n The results are the same as in [rLEC97c], Table 7!\n";

   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 4294967296, a = 1099087573. \n";
   maxdim = 8;
   m = 4294967296;
   m = m / 4;
   a.SetLength(2);
   a[1] = 1099087573;  
   FoMMRGLattice<NTL::ZZ, double>(m, maxdim, a, t);    
  std::cout << "\n The results are the same as in [rLEC97c], Table 1!\n";

  
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with lacunary indices\n";
   std::cout << "m = 2147483647, a = 16807, d =131072. \n";
   std::cout << "The lacunary indices are 1, 2, 3, d+1, d+2, d+3, 2*d+1, 2*d+2\n";  
   NTL::Vec<NTL::ZZ> lac; // Vector containing the lacunary indices.
   NTL::ZZ d;
   d = 131072;  
   lac.SetLength(8);
   lac[0] = 1;   lac[1] = 2;   lac[2] = 3;
   lac[3] = d+1;   lac[4] = d+2;   lac[5] = d+3;
   lac[6] = 2*d+1;   lac[7] = 2*d+2;
   m = 2147483647;
   a.SetLength(2);
   a[1] = 16807;
   t.SetLength(1);
   t[0] = 8;    
   FoMMRGLatticeLac<NTL::ZZ, double>(m, maxdim, a, t, lac);    
   std::cout << "\n The results are the same as in [rLEC97c], Table 2!\n";
}

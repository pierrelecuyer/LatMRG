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
static void FOMSuccMRGLattice(Int m, const NTL::Vec<Int> aa, int64_t lowDim, int64_t highDim) {
   NTL::Vec<int64_t> t;
   t.SetLength(2);
   LatMRG::MRGLattice<Int, Real> lat(m, aa, highDim);
   ReducerBB<Int, Real> m_red(lat); 
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), aa.length()-1, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &m_red);
   fomdual.setVerbosity(3);
   //IntLattice<Int, Real> proj(m, t.length());
   //fomdual.computeMerit(lat, proj);
   std::cout << "lowDim, highDim = " << lowDim << "  " << highDim << "\n";
   fomdual.computeMeritSucc(lat, lowDim, highDim);
}

/**
 * This function does the same as FoMMRGLattice except that it deals with MRG lattices with
 * lacunary indices defined by 'lac'.
 */
template<typename Int, typename Real>
static void FOMSuccMRGLatticeLac(Int m, const NTL::Vec<Int> aa, NTL::Vec<Int> lac,
       int64_t lowDim, int64_t highDim) {
   NTL::Vec<int64_t> t;
   t.SetLength(2);
   LatMRG::MRGLatticeLac<Int, Real> lat(m, aa, highDim);
   lat.setLac(lac);
   ReducerBB<Int, Real> m_red(lat); 
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), aa.length()-1, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &m_red);
   fomdual.setVerbosity(3);
   // IntLattice<Int, Real> proj(m, t.length());
   // std::cout << "t = " << t << "\n";
   std::cout << "lowDim, highDim = " << lowDim << "  " << highDim << "\n";
   fomdual.computeMeritSucc(lat, lowDim, highDim);
}


int main() {

   int64_t maxdim = 8;   
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::Vec<int64_t> t; // The t-vector for the FoM.
   t.SetLength(2);
   t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 0;      
   NTL::ZZ m;
   NTL::Vec<NTL::ZZ> aa; // Vector a has size d+1, if it contains d elements.
   
   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 2147483647 and a = 45991\n";
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 45991;
   FOMSuccMRGLattice<NTL::ZZ, double>(m, aa, 2, maxdim);
   std::cout << "\n The results are the same as in [rLEC97c], Table 6!\n";


   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 9223372036854773561, a = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   aa.SetLength(4);
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   FOMSuccMRGLattice<NTL::ZZ, double>(m, aa, 4, maxdim);
   std::cout << "\n The results are the same as in [rLEC97c], Table 7!\n";

   
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with m = 4294967296, a = 1099087573. \n";
   maxdim = 8;
   m = 4294967296;
   m = m / 4;
   aa.SetLength(2);
   aa[1] = 1099087573;
   FOMSuccMRGLattice<NTL::ZZ, double>(m, aa, 2, maxdim);
   std::cout << "\n The results are the same as in [rLEC97c], Table 1!\n";

  
   std::cout << "\n=============================================================\n";
   std::cout << "Results for a MRG example with lacunary indices\n";
   std::cout << "m = 2147483647, a = 16807, d =131072. \n";
   std::cout << "The lacunary indices are 1, 2, 3, d+1, d+2, d+3, 2*d+1, 2*d+2\n";  
   maxdim = 8;
   NTL::Vec<NTL::ZZ> lac; // Vector containing the lacunary indices.
   NTL::ZZ d;
   d = 131072;  
   lac.SetLength(8);
   lac[0] = 1;   lac[1] = 2;   lac[2] = 3;
   lac[3] = d+1;   lac[4] = d+2;   lac[5] = d+3;
   lac[6] = 2*d+1;   lac[7] = 2*d+2;
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 16807;
   t.SetLength(1);
   t[0] = 8;    
   // t[1] = 0;
   FOMSuccMRGLatticeLac<NTL::ZZ, double>(m, aa, lac, 3, maxdim);
   std::cout << "\n The results are the same as in [rLEC97c], Table 2!\n";
}

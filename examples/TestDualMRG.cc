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
#include "latticetester/ReducerStatic.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"

using namespace LatticeTester;

int main() {
    int64_t maxdim = 8;

    NTL::ZZ m;
    NTL::Vec<NTL::ZZ> a; // Vector a has size 4, a[j] contains a_j.

   a.SetLength(3);
   m = 101;
   a[1] = 37;
   a[2] = 22;

   LatMRG::MRGLattice<NTL::ZZ, double> lat(m, a, maxdim);

   lat.setUsePolynomialBasis(true);

   ReducerBB<NTL::ZZ, double> m_red(lat);   
   
   m_red.setVerbosity(2);

   lat.buildDualBasis(8);

   std::cout << " The dual basis is given by \n ############## \n" << lat.getDualBasis() << "\n";
   
}

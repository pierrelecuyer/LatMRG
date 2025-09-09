/**
 * This is a dummy program used to compile unused files of latmrg.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <algorithm>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGLatticeLac.h"
#include "latmrg/MLCGLattice.h"

using namespace LatticeTester;

int main() {

   std::cout << "\n=============================================================\n";
   std::cout << "TestCompile running.\n";

   NTL::Vec<int64_t> t;
   int64_t lowDim = 4;
   int64_t highDim = 15;
   int64_t k = 3;
   Int m(101);
   IntVec aa;
   aa.SetLength(4);
   aa[1] = 7;  aa[2] = 1;  aa[3] = 21;
   IntMat A;
   A.SetDims(3, 3);
   A[0][0] = 0;  A[0][1] = 1;  A[0][2] = 0;
   A[1][0] = 0;  A[1][1] = 0;  A[1][2] = 1;
   A[2][0] = 21;  A[2][1] = 1;  A[2][2] = 7;

   LatMRG::MRGLattice<Int, Real> mrg(m, aa, highDim);
   ReducerBB<Int, Real> red(mrg);
   WeightsUniform weights(1.0);
   NormaBestLat norma(log(m), k, highDim, L2NORM);
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red);
   fom.setVerbosity(3);
   std::cout << "\nUsing MRGLattice: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fom.computeMeritSucc(mrg, lowDim, highDim);

   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLattice, dual: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);

   A = (A * A) * A;
   LatMRG::MLCGLattice<Int, Real> mlcg(m, A, highDim, 3);
   red.setIntLattice(mlcg);
   std::cout << "\nUsing MLCGLattice: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   // mlcg.buildBasis(highDim);
   fom.computeMeritSucc(mlcg, lowDim, highDim);

   // NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   // FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MLCGLattice, dual: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mlcg, lowDim, highDim);

}

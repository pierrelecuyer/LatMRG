/**
 * This program illustrates various ways of using the MRG and MLCG lattice classes
 * to test the lattice structure of an MRG with or without lacunary indices.
 * In each case, the function `computeMeritSucc` computes the FOM for for successive
 * coordinates only, using the standard BKZ reduction followed by a BB for each
 * number of dimensions.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double
// With Int = int64_t, there are issues with computing characteristic polynomials.

#include <iostream>
#include <cstdint>
#include <algorithm>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latmrg/FlexModInt.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGLatticeLac.h"
#include "latmrg/MLCGLattice.h"
#include "latmrg/MLCGLatticeLac.h"

using namespace LatticeTester;

int main() {

   NTL::Vec<int64_t> t;
   int64_t lowDim = 4;
   int64_t highDim = 24;
   int64_t k = 3;
   Int m(101);

   IntVec aa;
   aa.SetLength(k+1);
   aa[1] = 7;  aa[2] = 1;  aa[3] = 21;

   IntMat A;
   A.SetDims(k, k);
   A[0][0] = 0;  A[0][1] = 1;  A[0][2] = 0;
   A[1][0] = 0;  A[1][1] = 0;  A[1][2] = 1;
   A[2][0] = 21;  A[2][1] = 1;  A[2][2] = 7;

   NTL::Mat<Int> Ak = A;
   for (int64_t i = 1; i < k; i++)
      Ak = Ak * A;
   for (int64_t i = 0; i < k; i++)
      for (int64_t j = 0; j < k; j++)
         Ak[i][j] = Ak[i][j] % m;

   std::cout << "\n=============================================================\n";
   std::cout << "TestMLCG101 running.\n";
   std::cout << "with k = " << k << ",  m = " << m << "\n";
   std::cout << "MRG with aa = " << aa << "\n";
   std::cout << "MLCG with A = \n" << A << "\n";
   std::cout << " and with Ak = A^k mod m = \n" << Ak << "\n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   std::cout << "\n=============================================================\n";

   IntVec lac;
   lac.SetLength(highDim);
   WeightsUniform weights(1.0);
   NormaBestLat norma(log(m), k, highDim, L2NORM);

   // Using MRG
   LatMRG::MRGLattice<Int, Real> mrg(m, aa, highDim);
   ReducerBB<Int, Real> red(mrg);
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red);
   fom.setVerbosity(3);
   std::cout << "\nUsing MRGLattice: \n";
   fom.computeMeritSucc(mrg, lowDim, highDim);

   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLattice, dual: \n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);


   // Using MLCG
   LatMRG::MLCGLattice<Int, Real> mlcg(m, Ak, highDim, k);
   red.setIntLattice(mlcg);
   std::cout << "\nUsing MLCGLattice with A^k: \n";
   mlcg.buildBasis(highDim);
   fom.computeMeritSucc(mlcg, lowDim, highDim);
   std::cout << "\nUsing MLCGLattice with A^k, dual: \n";
   fomdual.computeMeritSucc(mlcg, lowDim, highDim);


/*
   // Using MRG Lac
   LatMRG::MRGLatticeLac<Int, Real> mrglac(m, aa, highDim);
   for (int64_t i = 0; i < highDim; i++)
      lac[i] = i+1;
   mrglac.setLac(lac);
   red.setIntLattice(mrglac);
   std::cout << "\nUsing MRGLatticeLac, i_j=j: \n";
   mrglac.buildBasis(highDim);
   fom.computeMeritSucc(mrglac, lowDim, highDim);
   std::cout << "\nUsing MRGLatticeLac, dual: \n";
   fomdual.computeMeritSucc(mrglac, lowDim, highDim);

   // Using MLCG Lac with A^k
   LatMRG::MLCGLatticeLac<Int, Real> mlcglac3(m, Ak, highDim, 3);
   for (int64_t i = 0; i < highDim; i++)
      lac[i] = i + 1;
   mlcglac3.setLac(lac);
   red.setIntLattice(mlcglac3);
   std::cout << "\nUsing MLCGLatticeLac with Ak, i_j=j: \n";
   fom.computeMeritSucc(mlcglac3, lowDim, highDim);
   std::cout << "\nUsing MLCGLatticeLac dual: \n";
   fomdual.computeMeritSucc(mlcglac3, lowDim, highDim);
*/


   // Using MLCG Lac with A
   LatMRG::MLCGLatticeLac<Int, Real> mlcglac(m, k, highDim, k);
   for (int64_t i = 0; i < k; i++)
      lac[i] = i + 1;
   for (int64_t i = k; i < highDim; i++)
      lac[i] = k * (i - k + 2);
   mlcglac.setLac(lac);
   mlcglac.setA(A);
   red.setIntLattice(mlcglac);
   std::cout << "\nUsing MLCGLatticeLac with A: \n";
   fom.computeMeritSucc(mlcglac, lowDim, highDim);
   std::cout << "\nUsing MLCGLatticeLac dual: \n";
   fomdual.computeMeritSucc(mlcglac, lowDim, highDim);


   // Using MRG Again.
   std::cout << "\nUsing MRGLattice AGAIN \n";
   fom.computeMeritSucc(mrg, lowDim, highDim);
   std::cout << "\nUsing MRGLattice, dual: \n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);
}

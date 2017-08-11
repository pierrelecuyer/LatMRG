#include "latmrg/MMRGLattice.h"
#include "latmrg/PolyPE.h"

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/ntlwrap.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace NTL;
using namespace LatticeTester;

//============================================================================

namespace LatMRG
{

/* Max order for lacunary case in this class; takes too much memory.
 For order > ORDERMAX, use subclass MRGLatticeLac instead */
//PW_TODO à voir plus tard avec lacunary
#define ORDERMAX 100

//===========================================================================

MMRGLattice::MMRGLattice(const MScal & m, const MMat & A, int maxDim, int r,
                        NormType norm, LatticeType lat):
      IntLattice::IntLattice(m, r, maxDim, norm)
{
   m_A = A;
   m_latType = lat;
   m_lacunaryFlag = false;
   m_ip = new bool[1];
   init();
   //PW_TODO attention m_ip aussi initialisé dans init()

}


//===========================================================================

MMRGLattice::MMRGLattice(const MScal & m, const MMat & A, int maxDim, int r,
                       BVect & lac, NormType norm, LatticeType lat):

      IntLattice::IntLattice (m, r, maxDim, norm),
      //PW_TODO à décommenter ?
      //m_lac(lac, r),
      m_ip(0)
{
   m_A = A;
   m_latType = lat;
   m_lacunaryFlag = true;
   init();

   //initialization of B
   m_B.resize(lac.length(), A.size1());
   for (int k = 0; k < lac.length(); k++)
      m_B[k][conv<int>(lac[k]) - 1] = 1;
}


//===========================================================================

MMRGLattice::MMRGLattice(const MMRGLattice & lat):
  IntLattice::IntLattice (lat.m_modulo, lat.m_order,
                          lat.getDim(), lat.getNorm ()), m_lac(lat.m_lac)
{
   m_A = lat.m_A;
   m_latType = lat.m_latType;
   m_lacunaryFlag = lat.m_lacunaryFlag;

   m_ip = new bool[m_order];
   m_xi.SetLength (m_order);
   m_A.SetDims (m_order, m_order);
   m_sta.SetDims (m_order, m_order);

   int dim = getDim();
   int rmax = max(m_order, dim);
   m_wSI.SetDims (rmax, dim);

/*
   for (i = 0; i <= m_order; i++)
      m_xi[i] = lat.m_xi[i];
   for (i = 0; i <= m_order; i++)
      m_ip[i] = lat.m_ip[i];

   int j;
   for (i = 0; i <= m_order; i++)
      for (j = 0; j <= m_order; j++)
         m_sta[i][j] = lat.m_sta[i][j];

   for (i = 0; i <= maxDim; i++)
      for (j = 0; j <= maxDim; j++)
         m_wSI[i][j] = lat.m_wSI[i][j];
*/
   //PW_TODO : ça doit vraiment rester commenté ?
}

//============================================================================

MMRGLattice::~MMRGLattice ()
{
   kill();
}

//============================================================================

void MMRGLattice::kill()
{
   IntLattice::kill();
   if (0 != m_ip)
      delete[] m_ip;
   m_ip = 0;
   m_xi.kill();

   // PW_TODO : methode kill fonctionne sur une matrice ?
   //m_aCoef.kill();
   //m_A.kill();

   // PW_TODO à quoi ça sert de killer ça ?
   m_sta.kill();
   m_wSI.kill();
}

//=========================================================================

MMRGLattice & MMRGLattice::operator= (const MMRGLattice & lat)
{
   if (this == &lat)
      return *this;
   m_dim = lat.m_dim;
   copyBasis(lat);
   m_order = lat.m_order;
   m_ip = lat.m_ip;
   //m_shift = lat.m_shift;
   return *this;
   //MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
   //copy (lat);
   //return *this;
}

//===========================================================================

void MMRGLattice::init()
{
   kill(); //PW_TODO : wzf ?
   IntLattice::init();
   m_xi.SetLength(m_order);
   m_A.SetDims(m_order, m_order);
   if (m_order > ORDERMAX) {
      m_ip = new bool[1];
      m_sta.SetDims(1, 1);
   } else {
      m_ip = new bool[m_order];
      m_sta.SetDims(m_order, m_order);
   }
   int rmax = max(m_order, getDim());
   m_wSI.SetDims(rmax, getDim());

   double temp;
   conv(temp, m_modulo);
   double lgm2 = 2.0 * Lg (temp);
   calcLgVolDual2 (lgm2);
   //if (m_latType == ORBIT)
   //   initOrbit();
}

//=========================================================================

BScal & MMRGLattice::getLac (int j)
{
   if (isLacunary() && j <= m_lac.getSize() && j > 0)
      return m_lac.getLac(j);
   throw std::out_of_range("MMRGLattice::getLac");
}


//===========================================================================

void MMRGLattice::setLac(const Lacunary & lac)
{
   m_lac = lac;
   m_lacunaryFlag = true;
}


//===========================================================================

//PW_TODO merdasse à tester

string MMRGLattice::toStringGeneratorMatrix () const
{
   std::ostringstream out;
   out << "[";
   for (int i = 0; i < m_order; i++) {
      out << "[";
      for (int j = 0; j < (m_order-1); j++) {
         out << m_A[i][j] << " ";
      }
      if (i == (m_order-1))
         out << m_A[i][m_order-1] << "]";
      else
         out << m_A[i][m_order-1] << "]" << endl;
   }
   out << "]" << endl;

   return out.str ();
}


//===========================================================================

void MMRGLattice::buildBasis (int d)
{
   if (m_lacunaryFlag)
      buildLacunaryBasis(d, m_B);
   else
      buildNonLacunaryBasis(d);
}

//===========================================================================

void MMRGLattice::buildNonLacunaryBasis (int dimension)
// a basis is built in dimension d

{
   setDim(dimension);
   int sizeA = getOrder();
   m_basis.resize(dimension, dimension);
   m_vecNorm.resize(dimension);
   m_dualvecNorm.resize(dimension);
   setNegativeNorm ();
   setDualNegativeNorm ();

   // filling in the diagonal of m_basis
   for (int i = 0; i < sizeA; i++)
      m_basis[i][i] = 1;
   for (int i = sizeA; i < dimension; i++)
      m_basis[i][i] = m_modulo;

   // using genrator matrix A to complete the first lines of m_basis
   // with values generatred by the recurrence
   ZZ_p::init(m_modulo);
   MMatP temp;
   temp.SetDims(sizeA, sizeA);
   for (int i = 0; i < sizeA; i++)
      temp[i][i] = 1;

   int maxIter = floor(dimension/sizeA);

   for (int k = 1; k < maxIter+1; k++) {
      // calculation of transpose(A^k)
      temp *= conv<MMatP>(transpose(m_A));

      if (k == maxIter) { // we completed the end of m_basis matrix
         int residu = dimension - maxIter * sizeA;
         for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < residu; j ++)
               m_basis[i][k*sizeA +j] = conv<MScal>(temp[i][j]);
         }
      } else {
         for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < sizeA; j ++)
               m_basis[i][k*sizeA +j] = conv<MScal>(temp[i][j]);
         }
      }
   }

   // we create the dual lattice associated
   m_dualbasis.resize(dimension, dimension);
   CalcDual<BMat>(m_basis, m_dualbasis, dimension, m_modulo);

   if (!checkDuality())
      MyExit (1, "BUG in MMRGLattice::buildNonLacunaryBasis");
}


//===========================================================================

void MMRGLattice::buildLacunaryBasis (int dimension, BMat B)
{
   setDim(dimension);
   int sizeA = getOrder();
   m_basis.resize(dimension, dimension);
   m_vecNorm.resize(dimension);
   m_dualvecNorm.resize(dimension);
   setNegativeNorm ();
   setDualNegativeNorm ();

   // filling in the diagonal of m_basis
   for (int i = 0; i < sizeA; i++)
      m_basis[i][i] = 1;
   for (int i = sizeA; i < dimension; i++)
      m_basis[i][i] = m_modulo;

   // using genrator matrix A to complete the first lines of m_basis
   // with values generatred by the recurrence
   ZZ_p::init(m_modulo);
   MMatP temp;
   temp.SetDims(sizeA, sizeA);
   for (int i = 0; i < sizeA; i++)
      temp[i][i] = 1;

   int sizeB = B.NumRows();
   int maxIter = ceil( (dimension-sizeA) / sizeB );

   for (int k = 0; k < maxIter+1; k++) {
      // calculation of transpose(A^k)
      temp *= conv<MMatP>(transpose(m_A));


      if (k == maxIter) { // we completed the end of m_basis matrix
         int residu = dimension - sizeA - maxIter * sizeB;
         for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < residu; j ++)
               m_basis[i][sizeA + k*sizeB + j] = conv<MScal>( (temp * conv<MMatP>(transpose(B))) [i][j]);
         }
      } else {
         for (int i = 0; i < sizeA; i++) {
            for (int j = 0; j < sizeB; j ++)
               m_basis[i][sizeA + k*sizeB + j] = conv<MScal>( (temp * conv<MMatP>(transpose(B))) [i][j]);
         }
      }
   }

   // we create the dual lattice associated
   m_dualbasis.resize(dimension, dimension);
   CalcDual<BMat>(m_basis, m_dualbasis, dimension, m_modulo);

   if (!checkDuality())
      MyExit (1, "BUG in MMRGLattice::buildNonLacunaryBasis");

}

//===========================================================================

void MMRGLattice::getSubLine(MVect & vec, MMat& B, int lign, int jMin, int jMax)
{
    // both jMin and jMax are included
    vec.resize(jMax-jMin+1);
    for (int i = 0; i < (jMax-jMin+1); i++)
        vec[i] = B[lign][jMin+i];
}


//===========================================================================

void MMRGLattice::incDim()
{
   if (m_lacunaryFlag)
      incrementDimLacunaryBasis (m_B);
   else
      incrementDimNonLacunaryBasis ();
}

//===========================================================================

void MMRGLattice::incrementDimNonLacunaryBasis()
// X_n = A X_{n-1} mod m. We have: dimension >= order.
{
   IntLattice::incDim();
   int newDimension = getDim();
   int sizeA = getOrder();

   // ************* update of the primal lattice *************
   //  - we add a new coordinate to each vector v_i, this value being determined
   //    by the MMRG recurrence (even if the original vectors have been
   //    transformed linearly and we must apply the same transformations to their
   //    last coordinates).
   //  - we add an extra vector (0,..., 0, m) to complete this dimension
   //    increased basis.

   // we compute the number of steps required to reach the current state of
   // the generator for the considered dimension.
   int n = floor((newDimension-1) / sizeA);
   ZZ_p::init(m_modulo);
   MMatP temp;
   temp.SetDims(sizeA, sizeA);
   for (int i = 0; i < sizeA; i++)
      temp[i][i] = 1;
   for (int k = 0; k < n; k++)
     temp *= conv<MMatP>(transpose(m_A));
   // PW_TODO : could be useful to keep A^k in memory to shorten computation

   // update of the new v_i coordinates using the *temp* matrix and the first
   // coefficients of each line (can be seen as a seed vector). So this *temp*
   // matrix multiplied by this seed vector gives us the next values generated
   // by the MMRG for the considered dimension.
   MVect initialState;
   for (int i = 0; i < (newDimension-1); i++) {
     getSubLine(initialState, m_basis, i, 0, sizeA-1);
     initialState = conv<MVect>( transpose(temp) * conv<MVectP>(initialState) );
     m_basis[i][newDimension-1] = initialState[newDimension - n*sizeA -1];
   }
   m_basis[newDimension-1][newDimension-1] = m_modulo;


   // ************* update of the dual basis *************
   //  - we add a new 0 coordinate to each vector w_i in the dual basis.
   //  - for the new last line of the matrix, we add an extra vector,
   //    as described in L'Ecuyer's paper
   // PW_TODO : citer rédérence "Guide LatTester"

   MVect lastLine;
   lastLine.SetLength(newDimension);

   for (int i = 0; i < (newDimension-1); i++) {
      matrix_row<const BMat> row(m_dualbasis, i);
      lastLine -= m_basis[i][newDimension-1] * row;
   }
   for (int i = 0; i < (newDimension-1); i++)
      m_dualbasis[newDimension-1][i] = lastLine[i] / m_modulo;
   m_dualbasis[newDimension-1][newDimension-1] = 1;

   setNegativeNorm();
   setDualNegativeNorm();

   if (!checkDuality())
      MyExit (1, "BUG in MMRGLattice::incrementDimBasis");

}

//===========================================================================

void MMRGLattice::incrementDimLacunaryBasis(BMat B)
{
   IntLattice::incDim();
   int newDimension = getDim();
   int sizeA = getOrder();
   int sizeB = B.NumRows();

   // ************* update of the primal lattice *************
   //  - we add a new coordinate to each vector v_i, this value being determined
   //    by the MMRG recurrence (even if the original vectors have been
   //    transformed linearly and we must apply the same transformations to their
   //    last coordinates).
   //  - we add an extra vector (0,..., 0, m) to complete this dimension
   //    increased basis.

   // we compute the number of steps required to reach the current state of
   // the generator for the considered dimension.
   int n = ceil((newDimension-sizeA-1) / sizeB);
   ZZ_p::init(m_modulo);
   MMatP temp;
   temp.SetDims(sizeA, sizeA);
   for (int i = 0; i < sizeA; i++)
      temp[i][i] = 1;
   for (int k = 0; k < n+1; k++)
     temp *= conv<MMatP>(transpose(m_A));
   // PW_TODO : could be useful to keep A^k in memory to shorten computation

   // update of the new v_i coordinates using the *temp* matrix and the first
   // coefficients of each line (can be seen as a seed vector). So this *temp*
   // matrix multiplied by this seed vector gives us the next values generated
   // by the MMRG for the considered dimension.
   MVect initialState;
   MVect outputState;
   outputState.SetLength(sizeB);

   for (int i = 0; i < (newDimension-1); i++) {
      getSubLine(initialState, m_basis, i, 0, sizeA-1);
      outputState = conv<MVect>( conv<MMatP>(B) * transpose(temp) * conv<MVectP>(initialState) );
      m_basis[i][newDimension-1] = outputState[newDimension - sizeA - n*sizeB - 1];
   }
   m_basis[newDimension-1][newDimension-1] = m_modulo;


   // ************* update of the dual basis *************
   //  - we add a new 0 coordinate to each vector w_i in the dual basis.
   //  - for the new last line of the matrix, we add an extra vector,
   //    as described in L'Ecuyer's paper
   // PW_TODO : citer rédérence "Guide LatTester"

   MVect lastLine;
   lastLine.SetLength(newDimension);

   for (int i = 0; i < (newDimension-1); i++) {
      matrix_row<const BMat> row(m_dualbasis, i);
      lastLine -= m_basis[i][newDimension-1] * row;
   }
   for (int i = 0; i < (newDimension-1); i++)
      m_dualbasis[newDimension-1][i] = lastLine[i] / m_modulo;
   m_dualbasis[newDimension-1][newDimension-1] = 1;

   setNegativeNorm();
   setDualNegativeNorm();

   if (!checkDuality())
      MyExit (1, "BUG in MMRGLattice::incrementDimBasis");

}

//===========================================================================

} // end namespace

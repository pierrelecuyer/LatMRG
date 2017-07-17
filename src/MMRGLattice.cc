#include "latmrg/MMRGLattice.h"
#include "latmrg/PolyPE.h"

#include "latticetester/Util.h"
#include "latticetester/Const.h"

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
                       LatticeType lat, NormType norm):
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
                       BVect & lac, LatticeType lat, NormType norm):
      IntLattice::IntLattice (m, r, maxDim, norm), 
      m_lac(lac, maxDim), 
      m_ip(0)
{
   m_A = A;
   m_latType = lat;
   m_lacunaryFlag = true;
   init();
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
   m_A.kill();

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
   //kill(); //PW_TODO : wzf ?
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
      buildLacunaryBasis(d);
   else
      buildNonLacunaryBasis(d);
}

//===========================================================================

void MMRGLattice::buildNonLacunaryBasis (int dimension)
// a basis is built in dimension d

{
    int sizeA = m_A.NumCols();
    m_basis.resize(dimension, dimension);

    // filling in the diagonal of m_basis
    for (int i = 0; i < sizeA; i++)
        m_basis[i][i] = 1;
    for (int i = sizeA; i < dimension; i++)
        m_basis[i][i] = m_modulo;

   // using genrator matrix A to complete the first lines of m_basis
   // with values generatred by the recurrence
    ZZ_p::init(m_modulo);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++)
        temp[i][i] = 1;

    int maxIter = dimension/sizeA;

    for (int k = 1; k < maxIter+1; k++) {
        // calculation of transpose(A^k)
        temp *= conv<mat_ZZ_p>(transpose(m_A)); 

        if (k == maxIter) { // we completed the end of m_basis matrix
            int residu = dimension - maxIter * sizeA;
            for (int i = 0; i < sizeA; i++) {
                for (int j = 0; j < residu; j ++)
                    m_basis[i][k*sizeA +j] = conv<ZZ>(temp[i][j]);
            }
        } else {
            for (int i = 0; i < sizeA; i++) {
                for (int j = 0; j < sizeA; j ++)
                    m_basis[i][k*sizeA +j] = conv<ZZ>(temp[i][j]);
            }
        }
    }

    // we create the dual lattice associated
    m_dualbasis.resize(dimension, dimension);
    CalcDual<BMat>(m_basis, m_dualbasis, dimension, m_modulo);


    if (checkDuality ())
      cout << "*** Duality check: OK ***" << endl;

}


//===========================================================================

void MMRGLattice::buildLacunaryBasis (int d) 
{
   /*

   if (m_order > ORDERMAX)
      MyExit (1, "MRGLattice::buildLaBasis:   k > ORDERMAX");
   initStates();
   int IMax = m_lac.getSize();

   MVect b;
   b.SetLength(m_order);
   //Invert(m_aCoef, b, m_order);

   // b is the characteristic polynomial
   PolyPE::setM (m_modulo);
   PolyPE::setF(b);
   PolyPE pol;
   int ord = 0;

   // Construction d'un systeme generateur modulo m.
   for (int k = 0; k < IMax; k++) {
      // pour chaque indice lacunaire
      conv (m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod(m_e);
      pol.toVector (m_xi);

      ord = 0;
      for (int i = 1; i <= m_order; i++) {
         if (m_ip[i]) {
            ++ord;
            m_t5 = 0;
            for (int j = 1; j <= m_order; j++)
               m_t5 += m_sta[i][j] * m_xi[j - 1];
            m_wSI[ord][k] = m_t5;
         }
      }
   }

   //On veut s'assurer que la base m_v soit triangulaire (pour satisfaire les
   //conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur V_i >= i])
   //et de plein rang (on remplace les lignes = 0 par lignes avec m sur la
   //diagonale).
   Triangularization<BMat>(m_wSI, m_vSI, ord, IMax, m_modulo);
   CalcDual<BMat>(m_vSI, m_wSI, IMax, m_modulo);

   // Construire la base de dimension 1
   m_basis[0][0] = m_vSI[0][0];
   m_dualbasis[0][0] = m_wSI[0][0];
   setDim(1);

   setNegativeNorm();
   setDualNegativeNorm();

   for (int i = 2; i <= d; i++)
      incrementDimLacunaryBasis (IMax);

   */
}



//===========================================================================

void MMRGLattice::getSubLine(vec_ZZ & vec, mat_ZZ& B, int lign, int jMin, int jMax)
{
    // both jMin and jMax are included
    vec.SetLength(jMax-jMin+1);
    for (int i = 0; i < (jMax-jMin+1); i++)
        vec[i] = B[lign][jMin+i];
}


//===========================================================================

void MMRGLattice::incrementDim()
{
   if (m_lacunaryFlag)
      incrementDimLacunaryBasis (getDim());
   else
      incrementDimBasis ();
}

//===========================================================================

void MMRGLattice::incrementDimBasis()
// X_n = A X_{n-1} mod m. On a Dim >= Order.
{

    int oldDimension = m_basis.NumRows();
    int newDimension = oldDimension+1;
    int sizeA = m_A.NumRows();

    mat_ZZ primalBasisTemp = m_basis;
    m_basis.SetDims(newDimension, newDimension);

    for (int i = 0; i < oldDimension; i++) {
        for (int j = 0; j < oldDimension; j++)
            m_basis[i][j] = primalBasisTemp[i][j];
    }

    // calcul de la puissance a A en cours
    int n = floor(oldDimension / sizeA);
    ZZ_p::init(m_modulo);
    mat_ZZ_p temp;
    temp.SetDims(sizeA, sizeA);
    for (int i = 0; i < sizeA; i++) {
            temp[i][i] = 1;
    }

    // étape couteuse qui pourrait etre raccourcie en stockant A^k
    //--------------------------------------------------------------------
    for (int k = 1; k < n+1; k++) {
        temp *= conv<mat_ZZ_p>(transpose(m_A)); 
    }
    //--------------------------------------------------------------------

    // mise à jour de la nouvelle colonne de m_basis
    vec_ZZ initialState;
    for (int i = 0; i < oldDimension; i++) {
        getSubLine(initialState, m_basis, i, 0, sizeA-1);
        initialState = conv<vec_ZZ>( transpose(temp) * conv<vec_ZZ_p>(initialState) );
        m_basis[i][newDimension-1] = initialState[newDimension - n*sizeA -1];
    }

    m_basis[newDimension-1][newDimension-1] = m_modulo;

}

/*
   IntLattice::incDim();
   const int dim = getDim();
   //m_basis.setDim(dim);
   //m_w.setDim(dim);
   //write();

   for (int i = 0; i < dim; i++) {
      clear (m_vSI[0][i]);
      for (int j = 0; j < m_order; j++) {
         conv (m_t1, m_basis[i][dim - j - 2]);
         //m_t1 = m_t1 * m_aCoef[j];
         m_vSI[0][i] = m_vSI[0][i] + m_t1;
      }
      Modulo (m_vSI[0][i], m_modulo, m_vSI[0][i]);
      m_basis[i][dim-1] = m_vSI[0][i];
   }

   for (int i = 0; i < dim; i++)
      m_basis[dim-1][i] = 0;
   m_basis[dim-1][dim-1] = m_modulo;

   for (int i = 0; i < dim-1; i++)
      m_dualbasis[i][dim-1] = 0;
   m_dualbasis[dim-1][dim-1] = 1;

   for (int j = 0; j < dim-1; j++) {

      clear (m_t1);
      for (int i = 0; i < dim-1; i++) {
         m_t2 = m_dualbasis[i][j];
         m_t2 *= m_vSI[0][i];
         m_t1 -= m_t2;
      }
      Quotient (m_t1, m_modulo, m_t1);
      m_dualbasis[dim-1][j] = m_t1;
   }

   setNegativeNorm();
   setDualNegativeNorm();

   if (!checkDuality ())
      MyExit (1, "BUG");
}
*/










//===========================================================================

void MMRGLattice::incrementDimLacunaryBasis(int IMax)
{
   /*

   IntLattice::incDim();
   const int dim = getDim (); // new dimension (dim++)

   if (dim >= IMax) {
      MyExit (0,
    "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
   }

   for (int i = 0; i < dim-1; i++) {
      // v[i] -> VSI[0].
      for (int j = 0; j < dim-1; j++)
         m_vSI[0][j] = m_basis[i][j];
      clear (m_vSI[i][0]);

      for (int i1 = 0; i1 < dim-1; i1++) {
         ProdScal (m_vSI[0], m_wSI[i1], dim, m_wSI[i1][0]);
         Quotient (m_wSI[i1][0], m_modulo, m_wSI[i1][0]);
         m_t1 = m_wSI[i1][0] * m_vSI[i1][dim - 1];
         m_vSI[i][0] += m_t1;
      }
      Modulo (m_vSI[i][0], m_modulo, m_vSI[i][0]);
      m_basis[i][dim-1] = m_vSI[i][0];
   }

   for (int j = 0; j < dim-1; j++)
      m_basis[dim - 1][j] = 0;
   m_basis[dim -1][dim - 1] = m_vSI[dim -1][dim - 1];

   for (int i = 0; i < dim-1; i++)
      m_dualbasis[i][dim - 1] = 0;

   for (int j = 0; j < dim-1; j++) {
      clear (m_wSI[0][j]);
      for (int i = 0; i < dim-1; i++) {
         m_t1 = m_dualbasis[i][j];
         m_t1 *= m_vSI[i][0];
         m_wSI[0][j] += m_t1;
      }
      if (m_wSI[0][j] != 0)
         m_wSI[0][j] = -m_wSI[0][j];
      Quotient (m_wSI[0][j], m_vSI[dim - 1][dim - 1], m_wSI[0][j]);
      m_dualbasis[dim - 1][j] = m_wSI[0][j];
   }

   Quotient (m_modulo, m_vSI[dim - 1][dim - 1], m_t1);
   m_dualbasis[dim - 1][dim - 1] = m_t1;

   //setDim (dim + 1);
   setNegativeNorm ();
   setDualNegativeNorm ();

   */
}


//===========================================================================

}

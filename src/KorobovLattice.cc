#include "latticetester/Util.h"
#include "latmrg/KorobovLattice.h"
#include <cassert>

#ifdef WITH_NTL
using namespace NTL;
#else
using namespace boost::numeric::ublas;
#endif

using namespace std;
using namespace LatticeTester;


namespace LatMRG
{

KorobovLattice::KorobovLattice (const MScal & n, const MScal & a, int maxDim,
                                NormType norm) :
      LatMRG::IntLattice::IntLattice(n, 0, maxDim, norm)
{
   m_a = a;
   m_shift = 0;
   init();
}


//===========================================================================

KorobovLattice::KorobovLattice (const MScal & n, const MScal & a, int maxDim,
                                int t, NormType norm) :
      IntLattice::IntLattice(n, 0, maxDim, norm)
{
   m_a = a;
   m_shift = t;
   init();
}


//===========================================================================

KorobovLattice::~KorobovLattice()
{
}


//=========================================================================

KorobovLattice::KorobovLattice (const KorobovLattice & lat):
            IntLattice::IntLattice(lat.m_modulo, 0, lat.getDim (), lat.getNorm ())
{
   m_a = lat.m_a;
   m_shift = lat.m_shift;
}


//=========================================================================

KorobovLattice & KorobovLattice::operator= (const KorobovLattice & lat)
{
   if (this == &lat)
      return *this;
   m_dim = lat.m_dim;
   copyBasis(lat);
   m_order = lat.m_order;
   m_a = lat.m_a;
   m_shift = lat.m_shift;
   return *this;
}


//===========================================================================

void KorobovLattice::init()
{
//Erwan   IntLatticeBasis::init();
//   double temp;
//   conv (temp, m_m);
//   m_lgVolDual2[1] = 2.0 * Lg(temp);

//Erwan   for (int r = m_order + 1; r <= getMaxDim(); r++)
//Erwan      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}


//===========================================================================

std::string KorobovLattice::toStringCoef () const
{
   std::ostringstream out;
   out << m_a;
   return out.str ();
}


//===========================================================================

void KorobovLattice::buildBasis(int d)
{
   //assert(d <= getMaxDim());
   setDim(d);
   MScal tmp;
   conv(tmp, 1);

   for (int i = 0; i < m_shift; i++) {
      tmp *= m_a;
      tmp %= m_modulo;
   }

   for (int j = 1; j <= d; j++) {
      // V[1][j] = tmp % m;
      m_basis (0, j) = tmp;
      tmp *= m_a;
      tmp %= m_modulo;
   }

   for (int i = 2; i <= d; i++) {
      for (int j = 1; j <= d; j++) {
         if (i == j)
            m_basis (i, j) = m_modulo;
         else
            m_basis (i, j) = 0;
      }
   }

   // Build dual basis
   CalcDual<BMat>(m_basis, m_dualbasis, d, m_modulo);
}


//===========================================================================

void KorobovLattice::incDimSlow()
{
   // Temporaire: très lent. Reprogrammer.
   int d = getDim();
   buildBasis(d + 1);
   setNegativeNorm();
   setDualNegativeNorm();
}


//===========================================================================

void KorobovLattice::incDim()
{
   MScal tmp1, tmp2, tmp3; MVect vectmp1;// working variables
   IntLattice::incDim(); //Increment the dimenson of the lattice by 1
   const int dim = getDim(); //New dimension
   
   vectmp1.resize(dim);
   for (int i = 1; i < dim-1; i++) {
      conv (tmp2, m_basis (i, dim - 2));
      tmp1 = tmp2 * m_a;
      Modulo (tmp1, m_modulo, tmp1);
      m_basis (i, dim) = vectmp1(i) =  tmp1; //Erwan m_vSI (0, i) = tmp1;
   }

   matrix_row<BMat> row1(m_basis, dim - 2);
   SetZero (row1, dim - 2);

   for (int i = 0; i < dim-1; i++)
      m_basis (dim-1, i) = 0;
   m_basis (dim-1, dim-1) = m_modulo;

   for (int i = 0; i < dim-1; i++)
      m_dualbasis (i, dim-1) = 0;
   m_dualbasis (dim-1, dim-1) = 1;

   for (int j = 1; j < dim; j++) {
      clear (tmp1);
      for (int i = 1; i < dim; i++) {
         tmp2 = m_dualbasis (i, j);
         tmp2 *= vectmp1 (i);
         tmp1 -= tmp2;
      }
      Quotient (tmp1, m_modulo, tmp3);
      m_dualbasis (dim, j) = tmp3;
   }

   setNegativeNorm();
   setDualNegativeNorm();
}


//===========================================================================

}

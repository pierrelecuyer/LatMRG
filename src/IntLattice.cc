#include "latmrg/IntLattice.h"
#include "latticetester/Util.h"

#ifdef WITH_NTL
#include "NTL/quad_float.h"
#include "NTL/RR.h"
using namespace NTL;
#endif
using namespace std;
using namespace LatticeTester;


namespace LatMRG
{


//=========================================================================

IntLattice::IntLattice ( MScal modulo, int k, int maxDim, NormType norm ):
   IntLatticeBasis(maxDim, norm)
   {
      m_dim = maxDim;
      m_withDual = true;
      m_modulo = modulo;
      m_order = k;
      init ();
   }

//=========================================================================

IntLattice::IntLattice (const IntLattice & Lat):
   IntLatticeBasis(Lat)
   {
      m_order = Lat.m_order;
      init ();
   }

//=========================================================================


void IntLattice::init ()
{
   int dim = getDim ();
   kill ();
   double temp;
   conv (temp, m_modulo);
   m_lgVolDual2 = new double[dim];
   m_lgm2 = 2.0 * Lg (temp);
   m_lgVolDual2[0] = m_lgm2;
   m_vSI.resize(dim, dim);
   m_wSI.resize(dim, dim);
}

//=========================================================================

void IntLattice::kill ()
{
   if (m_lgVolDual2 == 0)
      return;
   delete [] m_lgVolDual2;
   m_lgVolDual2 = 0;
   m_vSI.clear();

   IntLatticeBasis::kill();
   if (!comp.empty()) {
      for (int s = 0; s < (int) comp.size(); s++)
         delete comp[s];
      comp.clear();
   }
}


//=========================================================================

IntLattice::~IntLattice ()
{
   kill ();
}

//=========================================================================

void IntLattice::incrementDimension ()
{
   IntLattice lattmp(*this);
   int dim = getDim();
   m_basis.resize(dim+1, dim+1);
   m_dualbasis.resize(dim+1, dim+1);
   m_vecNorm.resize(dim+1);
   m_dualvecNorm.resize(dim+1);

   for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
         m_basis(i,j) = lattmp.m_basis(i,j);
         m_dualbasis(i,j) = lattmp.m_dualbasis(i,j);
      }
      m_vecNorm(i) = lattmp.m_vecNorm(i);
      m_dualvecNorm(i) = lattmp.m_dualvecNorm(i);
   }
   setNegativeNorm(dim);
   setDualNegativeNorm(dim);
   setDim(dim+1);
}


void IntLattice::calcLgVolDual2 (double lgm2)
{
   int dim = getDim();
   int rmax = min(m_order, dim);

   m_lgVolDual2[0] = lgm2;
   for (int r = 1; r < rmax; r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1] + lgm2;
   // WARNING [David]: one version had `m_order` instead of `rmax`.
   // I am not sure which is the fix and which is the bug.
   for (int r = rmax; r < dim; r++)
      m_lgVolDual2[r] = m_lgVolDual2[r - 1];
}



}

#include "latmrg/IntLattice.h"

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
    }

//=========================================================================

IntLattice::IntLattice (const IntLattice & Lat):
    IntLatticeBasis(Lat)
    {
        m_order = Lat.m_order;
    }

//=========================================================================

void IntLattice::kill ()
{

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

void IntLattice::IncrementDimension ()
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



}

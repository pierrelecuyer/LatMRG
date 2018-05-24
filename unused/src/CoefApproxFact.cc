#include "latmrg/CoefApproxFact.h"
#include "latticetester/Util.h"

using namespace LatticeTester;


namespace LatMRG {

CoefApproxFact::CoefApproxFact(Zone *Z, MScal & m):
    Coefficient(), m_Z(Z), m_m(m)
{
}


//===========================================================================

CoefApproxFact::~CoefApproxFact()
{
   m_m = 0;
}


//===========================================================================

void CoefApproxFact::set (const MScal & q0, MVect & A, int & i)
{
   Zone::ZoneType No = m_Z->getNo ();
   MScal q = q0;
   if (q == 1)
      q = 2;
   if (m_Z->DivQ[No])
      Quotient (m_m, q, A[i]);
   else
      A[i] = q;
}

}

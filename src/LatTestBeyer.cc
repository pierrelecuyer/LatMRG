#include "latmrg/LatTestBeyer.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"
#include "latmrg/Merit.h"
#include "latticetester/Const.h"
#include "latticetester/Util.h"
#include "latmrg/PolyPE.h"
#include "latmrg/LatticeTest.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;


void testPrinting(BMat A, string name)
{
   for (int i = 0; i < A.size1(); i++) {
      for (int j = 0; j < A.size2(); j++)
         cout << name << "[" << i << "][" << j << "]=" << A(i,j) << "; ";
      cout << endl;
   }
}


namespace LatMRG
{

//===========================================================================

LatTestBeyer::LatTestBeyer (LatticeTester::IntLattice * lat): LatticeTest (lat)
{
   m_criter = BEYER;
}

//===========================================================================

bool LatTestBeyer::test (int fromDim, int toDim, double minVal[])
{
   m_fromDim = fromDim;
   m_toDim = toDim;
   init ();

   resetFromDim (m_lat->getOrder (), fromDim);
   while (m_lat->getDim () < fromDim)
      m_lat->incDim ();
   Reducer red (*m_lat);

   m_lat->dualize ();
   red.preRedDieter (0);
   m_lat->dualize ();

   while (true) {
      if (m_dualF)
         m_lat->dualize ();

      //cout << "CURRENT BASIS =\n";
      //testPrinting(m_lat->getBasis(), "testMatrix");

      bool success = red.reductMinkowski (0);
      int dim = m_lat->getDim ();
      if (success) {
         m_lat->updateScalL2Norm (0);
         m_lat->updateScalL2Norm (dim-1);

         double x1, x2;        // si VV[1] et VV[dim] sont tres
         // grands, il faudrait envisager de changer x1 et x2 en xdouble.
         conv (x1, m_lat->getVecNorm (0));
         conv (x2, m_lat->getVecNorm (dim-1));
         m_merit[dim] = x1 / x2;


         /*
         cout << "\nlac *** dim = " << dim << endl;
         cout << "lac *** x1 = " << x1 << endl;
         cout << "lac *** x2 = " << x2 << endl;
         cout << "lac *** m_merit = " << sqrt(m_merit[dim]) << endl;
         */


         // Si on sait deja que ce gen. ne pourra etre retenu,
         // on le rejette tout de suite et on arrete le test.
         if ((m_maxAllDimFlag && (m_merit[dim] < minVal[toDim]))
             || (m_merit[dim] < minVal[dim])) {
            m_merit[dim] = 0.0;
            return false;
         }
         if (3 == m_detailF) {
            dispatchLatUpdate(*m_lat);
         }

         prepAndDisp (dim);
         if (m_dualF)
            m_lat->dualize ();

      } else {
         m_merit[dim] = 0.0;
         return false;
      }

      if (dim == toDim)
         break;
      m_lat->incDim();
      red = Reducer(*m_lat);
   }

   return true;
}


//===========================================================================

void LatTestBeyer::init ()
{
   m_merit.setDim(m_toDim);
   const int N = 2;
   string header[N];
   header[0] = "q_t";
   header[1] = "Cumul CPU t(sec)";
   dispatchTestInit ("Beyer", header, N);
   timer.init();
}


//===========================================================================

void LatTestBeyer::prepAndDisp (int dim)
{
   const int N = 2;
   double results[N];
   results[0] = sqrt (m_merit[dim]);
   results[1] = timer.val (Chrono::SEC);
   dispatchResultUpdate (results, N);
}

//===========================================================================

}

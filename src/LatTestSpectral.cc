#include "latmrg/LatTestSpectral.h"
#include "latmrg/IntLattice.h"
#include "latticetester/Const.h"
#include "latmrg/Merit.h"
#include "latticetester/Util.h"
#include "latmrg/PolyPE.h"
#include "latmrg/LatticeTest.h"
#include "latticetester/Reducer.h"
#include <cmath>


//BOOST_DISPLAY
#include <boost/progress.hpp>

using namespace std;
using namespace NTL;
using namespace LatticeTester;

namespace LatMRG
{

//===========================================================================

LatTestSpectral::LatTestSpectral (Normalizer * normal,
              LatMRG::IntLattice * lat): LatticeTest (lat)
{
   m_criter = SPECTRAL;
   m_normalizer = normal;
}


//===========================================================================

LatTestSpectral::~LatTestSpectral ()
{
   m_boundL2.kill();
   delete [] m_S2toL2;
}


//===========================================================================

void LatTestSpectral::initLowerBoundL2 (int dim1, int dim2)
{
   if (m_lat->getNorm () == L2NORM) {
      for (int i = dim1; i <= dim2; i++) {
         m_S2toL2[i] = m_normalizer->getGamma (i);
         m_S2toL2[i] *= exp2 (m_lat->getLgVolDual2 (i) / i);
      }

   } else if (m_lat->getNorm () == L1NORM) {
      if (m_dualF) {
         for (int i = dim1; i <= dim2; i++)
            m_S2toL2[i] = trunc(exp2 ((m_lat->getLgVolDual2 (i) / 2.0
                          + m_normalizer->getGamma(i)) / i));
      } else {
         // Je ne suis pas sûr que ce soit correct pour le primal
         for (int i = dim1; i <= dim2; i++)
            m_S2toL2[i] = exp2 ((m_lat->getLgVolDual2 (i) / 2.0
                          + m_normalizer->getGamma(i)) / i);
      }

   } else {
      for (int i = dim1; i <= dim2; i++) {
         m_S2toL2[i] = 1.0;
      }
   }
}


//===========================================================================

void LatTestSpectral::setLowerBoundL2 (double S2)
{
   setLowerBoundL2(S2, 0);
}

void LatTestSpectral::setLowerBoundL2 (double S2, const double* weights)
{
   int i;
   NScal m2;
   conv(m2, m_lat->getModulo());
   m2 = m2 * m2;
   if (S2 <= 0.0) {
      for (i = m_fromDim; i <= m_toDim; i++)
         m_boundL2[i] = 0;
      return;
   }

   if (m_lat->getNorm () == L2NORM) {
      if (!m_dualF) {
         // On a multiplié les coordonnées du primal par m pour calculer
         // uniquement avec des entiers: meilleure précision.
         for (i = m_fromDim; i <= m_toDim; i++)
            conv(m_boundL2[i],m2 * S2 * m_S2toL2[i]);

      } else {
         for (i = m_fromDim; i <= m_toDim; i++)
            conv(m_boundL2[i], S2*m_S2toL2[i]);
      }

   } else {
      if (!m_dualF) {
         for (i = m_fromDim; i <= m_toDim; i++) {
            conv(m_boundL2[i], S2*m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            m_boundL2[i] = m_boundL2[i]*m_boundL2[i];
            // On a multiplié les coordonnées du primal par m pour calculer
            // uniquement avec des entiers: meilleure précision.
            m_boundL2[i] = m_boundL2[i] * m2;
         }

      } else {
         for (i = m_fromDim; i <= m_toDim; i++) {
            conv(m_boundL2[i], S2*m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            m_boundL2[i] = m_boundL2[i]*m_boundL2[i];
         }
      }
   }

   // Merit values will be scaled by the inverse weight after normalization.
   // It is the scaled values that need to be compared, but the Reducer doesn't know
   // about scaling. So we scale the bounds below by the weight, which ensures the
   // ordering of the scaled values is preserved.
   // FIXME: should it be the squared weight for L1NORM?
   if (weights) {
      for (i = m_fromDim; i <= m_toDim; i++)
         m_boundL2[i] *= weights[i];
   }
}


//===========================================================================

bool LatTestSpectral::test (int fromDim, int toDim, double minVal[])
{
   return test (fromDim, toDim, minVal, 0);
}

// weights is the array of the weights of all projections as follows:
//    weights[i] is the weight for projection [1, ..., i]
// weights == 0 means unit weight for all projections

//===========================================================================

bool LatTestSpectral::test (int fromDim, int toDim, double minVal[], const double* weights)
{

   //BOOST_DISPLAY
   boost::progress_display show_progress( toDim-fromDim+1 );

   m_merit.setDim(toDim);
   m_fromDim = fromDim;
   m_toDim = toDim;
   init ();
   resetFromDim (m_lat->getOrder (), fromDim);
   if (m_normalizer->getNorm () != m_lat->getNorm ()) {
      cout << "SpectralTest: conflict between NormType and Normalizer" <<
      endl;
      exit (EXIT_FAILURE);
   }

   double mr;
   conv (mr, m_lat->getModulo ());
   double temp;
   
   //NScal te;
   //double lgvv1 = 0.0;
   // "unused variables" according to the compiler
   
   while (m_lat->getDim () < fromDim)
      m_lat->incDim ();

   Reducer red (*m_lat);

   if (m_S2toL2[fromDim] <= 0.0)
   initLowerBoundL2 (fromDim, toDim);
   setLowerBoundL2 (minVal[toDim], weights);   // same S2 for all dim
   red.setBoundL2 (m_boundL2, fromDim, toDim);

   while (true) {

      //BOOST_DISPLAY
      ++show_progress;

      if (m_dualF)
         m_lat->dualize ();

      int dim = m_lat->getDim ();

      // pre-reduction step before BB with default parameters
      red.redBKZ(0.999999, 10, QUADRUPLE, dim);
      //PW_TODO: hard coding?

      //cout << "dual basis = \n" << m_lat->getBasis() << endl;
      //cout << "primal basis = \n" << m_lat->getDualBasis() << endl;

      if (red.shortestVector (m_lat->getNorm ())) {

         // Calcul de D2. Pour Norm # L2NORM, suppose que VV est a jour.
         if (m_lat->getNorm () == L2NORM) {
            m_lat->updateScalL2Norm (0);
            conv (temp, m_lat->getVecNorm (0));

            // to work with (0,1) variables otherwise we get the rescaled values
            /*if (!m_dualF) {
               conv(te, temp);
               NScal m2;
               conv(m2, m_lat->getModulo ());
               m2 = m2*m2;
               te = te / m2;
               conv(temp, te);
            }*/

         } else {
            conv (temp, red.getMinLength ());
            
            // to work with (0,1) variables otherwise we get the rescaled values
            /*if (!m_dualF)
               temp = temp / mr;
            */
         }
         if (3 == m_detailF) {
            dispatchLatUpdate(*m_lat);
            /*if (m_dualF) {
               dispatchLatUpdate (m_lat->getDualBasis());
               dispatchLatUpdate (m_lat->getBasis());
            } else {
               dispatchLatUpdate (m_lat->getBasis());
               dispatchLatUpdate (m_lat->getDualBasis());
            }*/
         }
         // weight factor:
         // require a higher merit for more important projections by dividing
         // their merit by the weight of the projection
         double weight = weights ? weights[dim] : 1.0;

         m_merit.getMerit (dim) = temp / weight;

         if (dim <= toDim) { // Calcul de S2.

            //cout << "dim = " << dim << endl;
            //cout << "order = " << m_lat->getOrder() << endl;
            //cout << "density AVANT = " << exp(m_normalizer->getLogDensity()) << endl;

            //updating value of matrix density.

      #if NTL_TYPES_CODE == 3
            if (dim <= m_lat->getOrder()) {
               if (m_dualF) // dual basis
                  m_normalizer->setLogDensity(conv<RScal>( - dim * log(m_lat->getModulo()) ));
               else // primal basis
                  m_normalizer->setLogDensity(conv<RScal>( dim * log(m_lat->getModulo())  ));
            }
      #else
            if (dim <= m_lat->getOrder()) {
               if (m_dualF) // dual basis
                  m_normalizer->setLogDensity( - dim * log(m_lat->getModulo()) );
               else // primal basis
                  m_normalizer->setLogDensity( dim * log(m_lat->getModulo()) );
            }
      #endif

            //cout << "density APRES = " << exp(m_normalizer->getLogDensity()) << endl;
            //cout << "m_normalizer->getBound = " << m_normalizer->getBound(dim) << endl;


            if (m_lat->getNorm () == L2NORM) {

               if (!m_dualF) { // in case we work with rescaled values
                  double normalizer = m_normalizer->getBound(dim);
                  m_merit[dim] = log(temp) - 2*log(normalizer) - 2*log(mr);
                  m_merit[dim] = exp(m_merit[dim]);
               } else { // general case
                  double normalizer = m_normalizer->getBound(dim);
                  m_merit[dim] = log(temp) - 2*log(normalizer);
                  m_merit[dim] = exp(m_merit[dim]);
               }

            } else if (m_lat->getNorm () == L1NORM) {

               if (!m_dualF) { // in case we work with rescaled values
                  double normalizer = m_normalizer->getBound(dim);
                  m_merit[dim] = log(temp) - log(normalizer) - log(mr);
                  m_merit[dim] = exp(m_merit[dim]);
               } else { // general case
                  double normalizer = m_normalizer->getBound(dim);
                  m_merit[dim] = log(temp) - log(normalizer);
                  m_merit[dim] = exp(m_merit[dim]);
               }

             } else
               m_merit[dim] = temp;

            m_merit[dim] /= weight;


            //cout << "m_merit[" << dim << "] = " << m_merit[dim] << endl;

            // Si on sait deja que ce gen. ne pourra etre retenu,
            // on le rejette tout de suite et on arrete le test.
            if ((m_maxAllDimFlag && m_merit[dim] < minVal[toDim])
                  || m_merit[dim] < minVal[dim]) {
           //    m_merit[dim] = 0.0;
               return false;
            }
         }
         else
            m_merit[dim] /= weight;

         prepAndDisp (dim);

      } else {
         m_merit[dim] = -1.0;
         return false;
      }

      if (dim == toDim)
         break;

      if (m_dualF)
            m_lat->dualize ();
      m_lat->incDim ();
      red = Reducer(*m_lat);
   }

   return true;
}


//===========================================================================

void LatTestSpectral::init ()
{
   m_merit.setDim(m_toDim);
   m_boundL2.SetLength (m_toDim+1);
   m_S2toL2 = new double[m_toDim+1];
   SetZero (m_S2toL2, m_toDim+1);
   const int N = 3;
   string header[N];
   if (m_lat->getNorm () == L2NORM) {
      if (m_invertF)
         header[0] = "d_t";
      else
         header[0] = "l_t";
      header[1] = "S_t";
      header[2] = "Cumul CPU t(sec)";
      dispatchTestInit ("SPECTRAL", header, N);

   } else {
      if (m_invertF)
         header[0] = "d_t";
      else
         header[0] = "N_t";
//      header[1] = "N_t^*";
      header[1] = "S_t";
      header[2] = "Cumul CPU t(sec)";
      dispatchTestInit ("SPECTRAL", header, N);
   }


   timer.init ();
}


//===========================================================================

void LatTestSpectral::prepAndDisp (int dim)
{
   const int N = 3;
   double results[N];

   if (2 == m_detailF)
      dispatchLatUpdate (*m_lat);
   else if (1 == m_detailF)
      dispatchLatUpdate (*m_lat, 0);

   if (m_lat->getNorm () == L2NORM) {
      results[0] = sqrt (m_merit.getMerit (dim));    // L_t
      if (m_invertF)
         results[0] = 1.0 / results[0];   // d_t = 1/L_t
      results[1] = sqrt (m_merit[dim]);

      results[2] = timer.val (Chrono::SEC);
      dispatchResultUpdate (results, N);

   } else {   // L1NORM
      results[0] = m_merit.getMerit (dim);
      if (m_invertF)
         results[0] = 1.0 / results[0];   // d_t = 1/N_t

#if 0
      // Calcule la valeur de N_t^*
      double x;
      conv (x, m_lat->getM());
      x = pow (x, (double) m_lat->getOrder());
      double y;
      if (m_dualF) {
         y = floor (pow (Factorial (dim) * x, 1.0 / dim));
      } else {
         x = 1.0/x;
         y = pow (Factorial (dim) * x, 1.0 / dim);
      }
      results[1] = y;
#endif

      results[1] = m_merit[dim];
      results[2] = timer.val (Chrono::SEC);
      dispatchResultUpdate (results, N);
   }
}

}

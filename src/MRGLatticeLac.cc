#include "latmrg/MRGLatticeLac.h"
#include "latmrg/PolyPE.h"
#include "latticetester/Util.h"

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

//===========================================================================

MRGLatticeLac::MRGLatticeLac (const MRGLatticeLac & lat): MRGLattice::
      MRGLattice (lat.m_modulo, lat.m_aCoef, lat.getDim (), lat.getOrder (),
                  lat.m_latType, lat.getNorm ()), m_lac (lat.m_lac)
{
   m_lacunaryFlag = true;
   MyExit (1, " MRGLatticeLac::  copy constructeur n'est pas terminé   ");
}


//=========================================================================

MRGLatticeLac & MRGLatticeLac::operator= (const MRGLatticeLac & lat)
{
   if (this == &lat)
      return * this;
   MyExit (1, " MRGLatticeLac::operator= n'est pas terminé   ");
   copyBasis (lat);
   return *this;
}


//===========================================================================

MRGLatticeLac::MRGLatticeLac (const MScal & m, const MVect & a, int maxDim,
                              int k, BVect & Lac, LatticeType lat,
                              NormType norm):
      MRGLattice::MRGLattice (m, a, maxDim, k, lat, norm),
      m_lac (Lac, maxDim)
{
   m_lacunaryFlag = true;
   m_sta.SetDims(1, 1);
   for (int i = 0; i < m_order; i++)
      m_aCoef[i] = a[i];
}


//============================================================================

MRGLatticeLac::~MRGLatticeLac ()
{
}


//=========================================================================

BScal & MRGLatticeLac::getLac (int j)
{
   if (j <= m_lac.getSize () && j > 0)
      return m_lac.getLac (j);
   throw std::out_of_range ("MRGLatticeLac::getLac");
}


//===========================================================================

void MRGLatticeLac::setLac (const Lacunary & lac)
{
   m_lac = lac;
   m_lacunaryFlag = true;
}


//===========================================================================

void MRGLatticeLac::buildBasis (int d)
{
   int ord = m_order;

   initStates ();
   int IMax = m_lac.getSize ();

   MVect b;
   b.SetLength (m_order+1);
   Invert (m_aCoef, b, m_order);

   // b is the characteristic polynomial
   PolyPE::setM (m_modulo);
   PolyPE::setF (b);
   PolyPE pol;

   // Construction d'un systeme generateur modulo m.
   for (int k = 0; k < IMax; k++) {
      // pour chaque indice lacunaire
      conv (m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod (m_e);
      pol.toVector (m_xi);

      for (int i = 0; i < m_order; i++) {
           m_wSI[i][k] = m_xi[i];
      }
   }

   //On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
   //les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
   //V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
   //avec m sur la diagonale).
    
   //Il serait possible de réserver m_wSI, m_vSI avec seulement IMax
   //lignes et colonnes quand order >> IMax. Mais il faudrait un nouveau
   //Triangularization qui pourrait nécessiter un appel pour chaque ligne,
   //mais qui sauverait beaucoup de mémoire.
   //Il n'est pas certain que cela en vaille la peine.
   Triangularization <BMat> (m_wSI, m_vSI, ord, IMax, m_modulo);
   CalcDual <BMat> (m_vSI, m_wSI, IMax, m_modulo);

   // Construire la base de dimension 1
   m_basis[0][0] = m_vSI[0][0];
   m_dualbasis[0][0] = m_wSI[0][0];
   setDim (1);

   setNegativeNorm ();
   setDualNegativeNorm ();

   for (int i = 1; i < d; i++)
      incDimBasis (IMax);

   // for debugging
   // trace("ESPION_2", i);
}


//===========================================================================

void MRGLatticeLac::incDimBasis (int IMax)
{
   MRGLattice::incDimLaBasis(IMax);
}


//===========================================================================

void MRGLatticeLac::initStates ()
/*
 * Initialise la matrice carrée Sta. La matrice Sta est d'ordre égal à l'ordre du
 * générateur. Elle contient un système de générateurs pour le groupe d'états
 * considérés.
 */
{
   clear (m_t2);

   if (m_latType == RECURRENT) {
      // check if a_k is relatively prime with m ==> m_t1 = 1
      m_t1 = GCD (m_aCoef[m_order], m_modulo);
      m_t1 = abs (m_t1);
      set9 (m_t2);
   }

   if (m_latType == FULL || m_latType == PRIMEPOWER || (m_t1 == m_t2)) {
      // m_sta is set to identity matrix
      calcLgVolDual2 (m_lgm2);

   } else {
      if (m_latType == ORBIT) {
         MyExit (1, "case ORBIT ne fonctionne pas");
      } else if (m_latType == RECURRENT) {
         MyExit (1, "case RECURRENT ne fonctionne pas");
      }
   }
}

//===========================================================================
}

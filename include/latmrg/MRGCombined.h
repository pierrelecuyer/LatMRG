#ifndef	LATMRG_MRGCOMBINED_H
#define	LATMRG_MRGCOMBINED_H

#include <NTL/mat_poly_ZZ.h>
#include <string>
#include "latticetester/FlexTypes.h"
#include "latmrg/EnumTypes.h"
// #include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

using namespace LatMRG;

/**
 * This class represents a combined MRG with two or more components of the same order \f$k\f$.
 * Each component \f$c\f$ has its own modulus \f$m_c\f$ and vector of multipliers \f$\mathbb{a}_c\f$,
 * and is represented internally as a `MRGComponent` object.
 * These components can be added or removed one by one as objects in a list.
 * The main tasks are to compute the modulus \f$m\f$, the vector of coefficients
 * \f$a_1,\dots,a_k\f$, and the period length of the combined MRG.
 */
template<typename Int> class MRGCombined {

public:

   /**
    * Constructor for a combined MRG with components of order `k`
    * given in a vector of `MRGComponent` objects.
    * They must be all of the same order `k`.
    */
   MRGCombined(int64_t k, std::vector<MRGComponent<Int>> &vcomp);

   /**
    * Constructs an `MRGCombined` object for components of order `k`.
    */
   MRGCombined(int64_t k);

   /**
    * Destructor.
    */
   ~MRGCombined();

   /**
    * Sets the vector of MRG components to `vcomp`.
    */
   void setComponents(std::vector<MRGComponent<Int>> &vcomp);

   /**
    * Adds component `comp` to the vector of components.
    */
   void addComponent(MRGComponent<Int> &comp);

   /**
    * Returns component `c`.
    */
   MRGComponent<Int> getComponent(int64_t c);

   /**
    * Returns the modulus of the combined MRG.
    */
   Int getModulus() const {
      return m_modulo;
   }

   /**
    * Returns the vector of multipliers for the combined MRG.
    */
   IntVec getaa() const {
      return m_aa;
   }

   /**
    * Returns the order \f$k\f$.
    */
   int64_t getOrder() const {
      return m_k;
   }

   /**
    * Computes the parameters of the combined MRG.
    */
   void computeCombination();


private:

   /**
    * The vector of MRG components.
    */
   std::vector<MRGComponent<Int>> m_vcomp;

   /**
    * The order \f$k\f$ of the recurrence.
    */
   int64_t m_k = 1;

   /**
    * The modulus \f$m\f$ of the combined MRG.
    */
   Int m_modulo;

   /**
    * The vector of multipliers \f$a_j\f$ of the combined MRG.
    */
   IntVec m_aa;

   /**
    * The vector of values of \f$n_c\f$ for the components.
    */
   IntVec m_nvec;

   /**
    * The period length \f$\rho\f$ for the combined MRG.
    * (Not used for now.)
    */
   Int m_rho;

};
// End class declaration

//===========================================================================
// IMPLEMENTTION

//=============================================================================
// Main constructor.
template<typename Int>
MRGCombined<Int>::MRGCombined(int64_t k) {
   m_k = k;
   m_aa.SetLength(k+1);
}

//===========================================================================
template<typename Int>
MRGCombined<Int>::~MRGCombined() {
   m_aa.kill();
}

//===========================================================================
template<typename Int>
void MRGCombined<Int>::setComponents(std::vector<MRGComponent<Int>> &vcomp) {
   m_vcomp = vcomp;
}

//===========================================================================
template<typename Int>
void MRGCombined<Int>::addComponent(MRGComponent<Int> &comp) {
   m_vcomp.push_back(comp);
}

//===========================================================================
template<typename Int>
MRGComponent<Int> MRGCombined<Int>::getComponent(int64_t c) {
   return m_vcomp[c];
}

//===========================================================================

template<typename Int>
void MRGCombined<Int>::computeCombination() {
   int64_t j, c;
   int64_t numcomp = m_vcomp.size();  // Number of components.
   m_nvec.SetLength(numcomp);
   NTL::Vec<Int> modvec;              // The individual moduli.
   modvec.SetLength(numcomp);
   NTL::Vec<Int> modvec1;             // Values of m/m_c.
   modvec1.SetLength(numcomp);
   m_modulo = Int(1);
   for (c = 0; c < numcomp; c++) {
      assert(m_vcomp[c].getOrder() == m_k);
      modvec[c] = m_vcomp[c].getModulus();
      m_modulo *= modvec[c];
   }
   // Compute the coefficients `a_j`, assuming that the m_c are relatively prime.
   for (j = 0; j < m_k; j++)  m_aa[j] = 0;
   for (c = 0; c < numcomp; c++) {
      modvec1[c] = m_modulo / modvec[c];
      // Beware: in PowerMod, the first arg must be smaller than the last one.
      m_nvec[c] = NTL::PowerMod(modvec1[c] % modvec[c], modvec[c] - 2, modvec[c]);
      // std::cout << " c = " << c << ",  m_nvec[c] = " << m_nvec[c] << "\n";
      Int prod = m_nvec[c] * modvec1[c];
      for (j = 1; j <= m_k; j++) {
         m_aa[j] += (m_vcomp[c].getaa())[j] * prod;
         m_aa[j] = m_aa[j] % m_modulo;
         // std::cout << " j = " << j << ",  m_aa[j] = " << m_aa[j] << "\n";
      }
   }
}

template class MRGCombined<std::int64_t> ;
template class MRGCombined<NTL::ZZ> ;

} // End namespace LatMRG
#endif

#ifndef	LATMRG_LCGCOMBINED_H
#define	LATMRG_LCGCOMBINED_H

#include <NTL/mat_poly_ZZ.h>
#include <string>
#include "latticetester/FlexTypes.h"
#include "latmrg/EnumTypes.h"
// #include "latmrg/IntFactorization.h"
// #include "latmrg/Primitivity.h"
#include "latmrg/LCGComponent.h"

namespace LatMRG {

using namespace LatMRG;

/**
 * This class represents a combined LCG with two or more components.
 * Each component \f$c\f$ has its own modulus \f$m_c\f$ and multiplier \f$\{a}_c\f$,
 * and is represented internally as a `LCGComponent` object.
 * These components can be added or removed one by one as objects in a list.
 * The main tasks are to compute the modulus \f$m\f$, the coefficient
 * \f$a\f$, and the period length of the combined LCG.
 */
template<typename Int> class LCGCombined {

public:

   /**
    * Constructor of a combined LCG.
    */
   LCGCombined();

   /**
    * Destructor.
    */
   ~LCGCombined();

   /**
    * Sets the vector of LCG components to `vcomp`.
    */
   // void setComponents(std::vector<LCGComponent<Int>> &vcomp);

   /**
    * Adds component `comp` to the vector of components.
    */
   void addComponent(LCGComponent<Int> &comp);

   /**
    * Returns component `c`.
    */
   LCGComponent<Int> getComponent(int64_t c);

   /**
    * Remove all the components.
    */
   void clearComponents();

   /**
    * Returns the modulus of the combined LCG.
    */
   Int getModulus() const {
      return m_modulo;
   }

   /**
    * Returns the multiplier of the combined LCG.
    */
   Int geta() const {
      return m_a;
   }

   /**
    * Computes the parameters of the combined LCG.
    */
   void computeCombination();


private:

   /**
    * The vector of LCG components.
    */
   std::vector<LCGComponent<Int>> m_vcomp;

   /**
    * The modulus \f$m\f$ of the combined LCG.
    */
   Int m_modulo;

   /**
    * The multiplier `a` of the combined LCG.
    */
   Int m_a;

   /**
    * The vector of values of \f$n_c\f$ for the components,
    * after computing the parameters of the combination.
    */
   IntVec m_nvec;

   /**
    * The period length \f$\rho\f$ for the combined LCG.
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
LCGCombined<Int>::LCGCombined() {}

//===========================================================================
template<typename Int>
LCGCombined<Int>::~LCGCombined() {}

//===========================================================================
template<typename Int>
void LCGCombined<Int>::addComponent(LCGComponent<Int> &comp) {
   m_vcomp.push_back(comp);
}

//===========================================================================
template<typename Int>
LCGComponent<Int> LCGCombined<Int>::getComponent(int64_t c) {
   return m_vcomp[c];
}

//===========================================================================
template<typename Int>
void LCGCombined<Int>::clearComponents() {
   m_vcomp.clear();
   // clear(m_nvec);
}

//===========================================================================

template<typename Int>
void LCGCombined<Int>::computeCombination() {
   int64_t c;
   int64_t numcomp = m_vcomp.size();  // Number of components.
   m_nvec.SetLength(numcomp);
   NTL::Vec<Int> modvec;              // The individual moduli.
   modvec.SetLength(numcomp);
   NTL::Vec<Int> modvec1;             // Values of m/m_c.
   modvec1.SetLength(numcomp);
   m_modulo = Int(1);
   for (c = 0; c < numcomp; c++) {
      modvec[c] = m_vcomp[c].getModulus();
      m_modulo *= modvec[c];
   }
   // Compute the coefficients `a_j`, assuming that the m_c are relatively prime.
   m_a = Int(0);
   for (c = 0; c < numcomp; c++) {
      modvec1[c] = m_modulo / modvec[c];
      // Beware: in PowerMod, the first arg must be smaller than the last one.
      m_nvec[c] = NTL::PowerMod(modvec1[c] % modvec[c], modvec[c] - 2, modvec[c]);
      // std::cout << " c = " << c << ",  m_nvec[c] = " << m_nvec[c] << "\n";
      // Int prod = m_nvec[c] * modvec1[c];
      m_a += m_vcomp[c].geta() * m_nvec[c] * modvec1[c];
      m_a = m_a % m_modulo;
   }
}

template class LCGCombined<std::int64_t> ;
template class LCGCombined<NTL::ZZ> ;

} // End namespace LatMRG
#endif

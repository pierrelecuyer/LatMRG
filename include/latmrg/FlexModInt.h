#ifndef LATMRG_FLEXMODINT_H
#define LATMRG_FLEXMODINT_H

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
//#include "NTL/ZZ_pXFactoring.h"
#include "NTL/ZZ_pE.h"
#include "NTL/lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
//#include "NTL/lzz_pXFactoring.h"

namespace LatMRG {

// The default case: use NTL::ZZ
template<typename Int>
class ModInt {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};

// Specialization for the case int64_t
template<>
class ModInt<std::int64_t> {
public:
   typedef NTL::zz_p IntP;
   typedef NTL::zz_pE PolE;
   typedef NTL::vec_zz_p IntVecP;
   typedef NTL::mat_zz_p IntMatP;
   typedef NTL::zz_pX PolX;
};
/*
template<>
class ModInt<NTL::ZZ> {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};
*/

/**
 * This file provides
 * The type `Int` is used to represent the polynomial coefficients.
 * The possible associated types `IntVec` are `int64_t*` and <tt>vec_ZZ</tt>.
 * The type `PolE` for the polynomials may be chosen
 * as <tt>zz_pE</tt> when \f$m\f$ is small enough, or may be set to
 * <tt>ZZ_pE</tt> which is implemented with the big integer type <tt>NTL::ZZ_p</tt>.
 */


/**
 * Sets to `m` the modulus used by NTL for its `IntP` calculations.
 */
//template<typename Int>
//static void setModulus(const Int &m);

//===========================================================================
// Implementation

template<typename Int>
static void setModulus(const Int &m) {
   ModInt<Int>::IntP::init(m);
}

}// namespace LatMRG
#endif

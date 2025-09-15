#ifndef LATMRG_FLEXMODINT_H
#define LATMRG_FLEXMODINT_H

//#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"  // Integers modulo p
#include "NTL/ZZ_pX.h" // Polynomials with coefficients modulo p
#include "NTL/ZZ_pE.h" // Polynomials modulo a polynomial P(z), modulo p.
//#include "NTL/ZZ_pXFactoring.h"
#include "NTL/lzz_p.h"
#include "NTL/lzz_pX.h"
#include "NTL/lzz_pE.h"
//#include "NTL/lzz_pXFactoring.h"

namespace LatMRG {

/**
 * The `FlexModInt` class permits one to use the NTL functions that work on
 * integers, vectors, and polynomials over `Z_p` (with arithmetic modulo p)
 * using the flexible type `Int` for the modulus and coefficients.
 * In NTL, these functions are in different files, with names `ZZ_*` for ZZ
 * and `zz_*` for `int64_t`. Below, we give the generic (template) names
 * `IntP, PolE, PolX, IntVecP, IntMatP` to some of these files, and the meaning
 * of these names depend on the template `Int`.
 * To use the appropriate version for type `Int`, it suffices to prefix the
 * generic file name by `FlexModInt<Int>::`.  For example,
 * `FlexModInt<Int>::IntP::init(p)` calls the `init` function from `ZZ_p` or `zz_p`,
 * and `FlexModInt<Int>::PolX f;` declares a polynomial from `ZZ_pX` or `zz_pX`,
 * depending on the meaning of `Int`.
 * The modulus `p` that is used in those classes for a given `Int` type must be
 * set via `setModulusIntP<Int>(p);`.
 */


// The default case: use NTL::ZZ
template<typename Int>
class FlexModInt {
public:
   typedef NTL::ZZ_p IntP;   // Integers modulo p
   typedef NTL::ZZ_pX PolX;  // Polynomials with coeff. in Z_p
   typedef NTL::ZZ_pE PolE;  // Polynomials modulo P(z), with coeff. in Z_p
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   static NTL::ZZ_p to_Int_p(Int a) {return NTL::to_ZZ_p(a);};
   static void mod_init(Int m) {NTL::ZZ_p::init(m);}
};

// Specialization for the case int64_t
template<>
class FlexModInt<std::int64_t> {
public:
   typedef NTL::zz_p IntP;
   typedef NTL::zz_pX PolX;
   typedef NTL::zz_pE PolE;
   typedef NTL::vec_zz_p IntVecP;
   typedef NTL::mat_zz_p IntMatP;
   static NTL::zz_p to_Int_p(int64_t a) {return NTL::to_zz_p(a);};
   static void mod_init(int64_t m) {NTL::zz_p::init(m);}
};

/*
// Not needed, I think.
template<>
class FlexModInt<NTL::ZZ> {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pX PolX;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   static NTL::ZZ_p to_Int_p(NTL::ZZ a) {return NTL::to_ZZ_p(a);};
   static void mod_init(NTL::ZZ m) {NTL::ZZ_p::init(m);}
};
*/

/**
 * Sets to `m` the modulus p used by NTL for its `IntP` calculations.
 */
template<typename Int>
static void setModulusIntP(const Int &m);

//===========================================================================
// Implementation

/**
 * Sets to `p` the modulus used by NTL for its `IntP` calculations.
 */
template<typename Int>
static void setModulusIntP(const Int &p) {
   FlexModInt<Int>::IntP::init(p);
}

} // namespace LatMRG
#endif

#ifndef	LATMRG_MRGCOMPONENT_H
#define	LATMRG_MRGCOMPONENT_H

// #include "latmrg/Modulus.h"
#include <NTL/mat_poly_ZZ.h>
#include <string>
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Primitivity.h"

namespace LatMRG {

using namespace LatMRG;

/**
 * This class offers tools to verify if an MRG recurrence or order \f$k\f$ modulo \f$m\f$,
 * of the form
 * \f[
 *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mbox{ mod } m,
 * \f]
 * has full period or not.  The object is constructed for given values of
 * \f$k\f$ and \f$m\f$, and can be used for several values of the
 * coefficients \f$a_1,\dots,a_k\f$.
 * It stores information such as the factorizations of
 * \f$m-1\f$ and of \f$r = (m^k-1)/(m-1)\f$, which are needed to check the full period
 * conditions.  It does not look at the lattice structure.
 */
template<typename Int> class MRGComponent {

public:

   /**
    * Constructor for modulus \f$m\f$ and order \f$k\f$.
    * The vector of multipliers must be set separately by `setaa`.
    * The Arguments `decom1` and `decor` specify the type of the prime factor
    * decomposition of \f$m-1\f$ and \f$r=(m^k-1)/(m-1)\f$, respectively.
    * The type `DecompType` is defined in `EnumTypes` and offers the
    * following choices:  `DECOMP`: the integer will be factored,
    * `DECOMP_WRITE`: it will be factored and the prime factors written in a file,
    * `DECOMP_READ`: the prime factors are read from a file,
    * `DECOMP_PRIME`: the integer is assumed to be prime.
    * The names `filem1` and `filer` (character strings) specify the file where the factors of
    * \f$m-1\f$ and \f$r\f$, are read or written when the type
    * <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt> is used.
    * The given files must be accessible by the program.
    * The format for the factorization in these files is described in class `IntFactorization`.
    */
   MRGComponent(const Int &m, int64_t k, DecompType decompm1, const char *filem1, DecompType decompr,
         const char *filer);

   /**
    * Same as the previous constructor, but with \f$m = b^e + c\f$.
    */
   MRGComponent(Int b, int64_t e, Int c, int64_t k, DecompType decompm1, const char *filem1,
         DecompType decompr, const char *filer);

   /**
    * Destructor.
    */
   ~MRGComponent();

   /**
    * Sets the vector of multipliers to `aa`, with `aa[j]` containing \f$a_j\f$.
    * The order \f$k\f$ of the MRG is set equal to the length of this vector, plus 1.
    * If we only want to call `maxPeriod`, there is no need to call `setaa` before.
    */
   void setaa(const IntVec &aa);

   /**
    * Returns the vector of multipliers of the recurrence, in same format as `aa` above.
    */
   IntVec getaa() const {
      return m_aa;
   }

   /**
    * Returns `true` if the coefficients \f$aa\f$ give an MRG with maximal
    * period; returns `false` otherwise.
    */
   bool maxPeriod(const IntVec &aa);

   /**
    * Assumes that the maximal-period condition (1) holds and checks only for
    * conditions 2 and 3. Returns `true` iff these two conditions hold for the
    * vector of coefficients \f$aa\f$.
    */
   bool maxPeriod23(const IntVec &aa);

   /**
    * Returns the value of the modulus \f$m\f$ of the recurrence.
    */
   const Int getModulus() const {
      return m_m;
   }

   /**
    * Returns the value of the order \f$k\f$ of the recurrence.
    */
   int64_t getOrder() const {
      return m_k;
   }

   /**
    * Returns the value of `b` for which `b^e+r = m`, when available.
    */
   Int& getb() {
      return m_b;
   }

   /**
    * Returns the value of `e` for which `b^e+r = m`, when available.
    */
   int64_t& gete() {
      return m_e;
   }

   /**
    * Returns the value of `c` for which `b^e+c = m`, when available.
    */
   Int& getc() {
      return m_c;
   }

   /**
    * Returns a const reference to the initial state `orbitSeed`, when available.
    * */
   IntVec& getOrbitSeed() {
      return orbitSeed;
   }

   /**
    * Returns a descriptor of this object as a string.
    */
   std::string toString();

   /**
    * Sets the type of this generator component.    ??????
    void setType(GenType type) {
    m_type = type;
    }
    * Gets the type of this component.
    *
    GenType getType() {
    return m_type;
    }
    */

private:

   /**
    * This is called by the constructor, with the same arguments.
    */
   void init(const Int &m, int64_t k, DecompType decom1, const char *filem1, DecompType decor,
         const char *filer);

   /**
    * The type of generator this stores. Should be MRG, MMRG or MWC. The
    * default value is LCG.
    * */
   // GenType m_type = LCG;   // This should not be here, but should depend on the subclass, it seems.
   /**
    * The prime factor decomposition of \f$m-1\f$.
    */
   IntFactorization<Int> ifm1;

   /**
    * The prime factor decomposition of \f$m\f$.
    * Used to compute the full period of a LCG with a carry.   // Not here!    *******
    */
   // IntFactorization<Int> factor;
   /**
    * The prime factor decomposition of \f$m-1\f$.
    * This is used when computing the full period of a matrix generator.
    * This is because we use NTL to compute the characteristic polynomial
    * of the matrix since it is already implemented.
    */
   // IntFactorization<NTL::ZZ> ifm2;
   /**
    * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
    * \f$k\f$ is the order of the recurrence.
    */
   IntFactorization<Int> ifr;

   /**
    * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
    * \f$k\f$ is the order of the recurrence.
    * This is used when computing the full period of a matrix generator.
    * This is because we use NTL to compute the characteristic polynomial
    * of the matrix since it is already implemented.
    */
   // IntFactorization<NTL::ZZ> ifr2;
   /**
    * The modulus \f$m\f$ of the recurrence, as a  `Modulus` object.
    */
   // Modulus<Int> modulus;    //  Want this ????
   /**
    * The modulus \f$m\f$ of this MRG.
    */
   Int m_m;

   /**
    * The order \f$k\f$ of the recurrence.
    */
   int64_t m_k = 1;

   /**
    * The multipliers \f$a_i\f$ of the recurrence, \f$i = 1, ..., k\f$.
    */
   IntVec m_aa;

   /**
    * The generator matrix \f$A\f$ of the recurrence for MMRG   ******  ???
    */
   // IntMat m_A;    // Want it here?
   /**
    * Basis `b` with `b^e+c = m`.
    * */
   Int m_b;

   /**
    * Exponent `e` with `b^e+c = m`.
    * */
   int64_t m_e;

   /**
    * Rest `c` with `b^e+c = m`.
    * */
   Int m_c;

   /**
    * The period length \f$\rho\f$ for this MRG. For now, the
    * recurrence is assumed to correspond to a primitive polynomial and
    * `rho` is calculated as
    * \f[
    *   \rho= m^k - 1
    * \f]
    * This value is calculated by `MRGLatticeFactory` and stored here for
    * simplicity.                     Used?   *****************
    */
   Int rho;

   /**
    * Value needed for the calculation of the multipliers of a combined
    * MRG. It is defined by
    * \f[
    *   n_j = (m/m_j)^{-1} \mbox{ mod } m_j \qquad\mbox{ for } j = 1,…,J,
    * \f]
    * where \f$n_j = \f$ `nj`, \f$m_j\f$ is this object’s modulus `m`,
    * \f$m\f$ is the calculated modulus for the combined MRG (see class
    * <tt>MRGLatticeFactory</tt>), and \f$(m/m_j)^{-1} \mbox{ mod } m_j\f$
    * is the inverse of \f$m/m_j\f$ modulo \f$m_j\f$. This value is
    * calculated by `MRGLatticeFactory` and stored here for simplicity.
    */
   Int nj;     // Not sure if we need this here.   ***************

   /**
    * Contains the initial state of the cycle we want to analyze when the
    * lattice type is `ORBIT`. It is made of \f$k\f$ numbers.
    *                            *********   But this class does not look at lattices!
    */
   IntVec orbitSeed;

};
// End class declaration

//===========================================================================
// IMPLEMENTTION

=============================================================================
// Main constructor.
template<typename Int>
MRGComponent<Int>::MRGComponent(const Int &m, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   m_b = m;
   m_e = 1;
   m_c = Int(0);
   init(m, k, decompm1, filem1, decompr, filer);
}

//===========================================================================

// Constructor with alternative format for m.
template<typename Int>
MRGComponent<Int>::MRGComponent(Int b, int64_t e, Int c, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   Int m = NTL::power(b, e) + c;
   m_b = b;
   m_e = e;
   m_c = c;
   init(m, k, decompm1, filem1, decompr, filer);
}

//===========================================================================

/*
 // Copy constructor.
 template<typename Int>
 MRGComponent<Int>::MRGComponent(const MRGComponent<Int> &lat) :
 ifm1(lat.ifm1), ifr(lat.ifr), m_k(lat.m_k) {
 m_m = lat.m_m;
 nj = lat.nj;
 rho = lat.rho;
 //   a.kill();
 m_aa.resize(m_k);
 m_aa = lat.m_aa;
 //   orbitSeed.kill();
 orbitSeed.resize(m_k);
 orbitSeed = lat.orbitSeed;

 // m_type = lat.m_type;
 //if (m_type == MWC)
 //    m_MWCb = lat.m_MWCb;
 m_b = lat.m_b;
 m_e = lat.m_e;
 m_c = lat.m_c;
 }

 //===========================================================================

 // Assignment constructor.
 template<typename Int>
 MRGComponent<Int>& MRGComponent<Int>::operator=(const MRGComponent<Int> &lat) {
 if (this != &lat) {
 m_k = lat.m_k;
 m_m = lat.m_m;
 nj = lat.nj;
 rho = lat.rho;
 //    a.kill();
 m_aa.resize(m_k);
 m_aa = lat.m_aa;
 //     orbitSeed.kill();
 orbitSeed.resize(m_k);
 orbitSeed = lat.orbitSeed;
 //      ifm1 = lat.ifm1;
 //      ifr = lat.ifr;
 }
 // m_type = lat.m_type;
 return *this;
 }
 */

//============================================================================
template<typename Int>
void MRGComponent<Int>::init(const Int &m, int64_t k, DecompType decompm1, const char *filem1,
      DecompType decompr, const char *filer) {
   m_m = m;
   setModulusIntP<Int>(m_m);
   m_k = k;
   m_aa.SetLength(m_k);
   orbitSeed.SetLength(m_k);   // Used ???

   Int mm1;  // m-1
   mm1 = m - 1;
   ifm1 = IntFactorization<Int>(mm1);
   ifm1.decompToFactorsInv (decompm1, filem1);
   Int r;
   r = (NTL::power(m, m_k) - 1) / (m - 1);
   ifr = IntFactorization<Int>(r);
   ifr.decompToFactorsInv (decompr, filer);
}

//===========================================================================

/**
 template<typename Int>
 MRGComponent<Int>::MRGComponent(Modulus<Int> &modu, int64_t k, DecompType decom1,
 const char *filem1, DecompType decor, const char *filer) {
 init(modu.m, k, decom1, filem1, decor, filer);

 LatticeTester::PrimeType status = LatticeTester::IntFactor < Int
 > ::isPrime(modu.m, 100);
 if (status == LatticeTester::PRIME || LatticeTester::PROB_PRIME == status) {
 modu.primeF = true;
 } else {
 std::cout << " WARNING:  m is NOT prime" << std::endl;
 modu.primeF = false;
 }
 module = modu;
 }
 */

//===========================================================================
template<typename Int>
MRGComponent<Int>::~MRGComponent() {
   //a.kill();
   //orbitSeed.kill();
}

//===========================================================================

template<typename Int>
void MRGComponent<Int>::setaa(const IntVec &aa) {
   // m_aa.SetLength(aa.length());
   m_aa = aa;
}

//===========================================================================

/*
 template<typename Int>
 void MRGComponent<Int>::setA(const IntMat &A) {
 m_A.SetDims(A.NumRows(), A.NumCols());
 m_A = A;
 }
 */

//===========================================================================

template<typename Int>
bool MRGComponent<Int>::maxPeriod(const IntVec &aa) {
   return isPrimitive(aa, m_m, ifm1, ifr);
}

//===========================================================================
// This is for a matrix LCG, does not belong here.   *************
/*

 template<typename Int>
 bool MRGComponent<Int>::maxPeriod(const IntMat &a0) {
 // Converting everything to NTL::ZZ to ease the characteristic polynomial
 // computation
 Primitivity<NTL::ZZ>::setM(NTL::ZZ(getM()));
 NTL::ZZX poly;
 NTL::matrix<NTL::ZZ> mat(m_k, m_k);
 for (int64_t i = 0; i < m_k; i++) {
 for (int64_t j = 0; j < m_k; j++) {
 mat[i][j] = NTL::ZZ(a0[i][j]);
 }
 }
 // Characteristic polynomial computation
 NTL::CharPoly(poly, mat);
 // Copying the polynomial to a vector
 NTL::vector<NTL::ZZ> vec(NTL::VectorCopy(poly, m_k + 1));
 Primitivity<NTL::ZZ>::setF(vec);
 Primitivity<NTL::ZZ> pol;
 PrimitiveInt < NTL::ZZ > privfm(ifm2, NTL::ZZ(getM()), 1);
 return pol.isPrimitive(privfm, ifr2);
 }
 */

//===========================================================================
// Must be rewritten!
template<typename Int>
bool MRGComponent<Int>::maxPeriod23(const IntVec &aa) {
   return isPrimitive(aa, m_m, ifm1, ifr);
}

//===========================================================================

template<typename Int>
std::string MRGComponent<Int>::toString() {
   std::ostringstream os;
   os << "MRGComponent:";
   os << "\n   m = " << m_m << " = " << m_b << "^" << m_e << " + " << m_c;
   os << "\n   k = " << m_k;
   os << "\n   a = ";
   std::string str(os.str());
   std::string s2 = LatticeTester::toString(m_aa, m_k);
   str += s2;
   str += "\n";
   return str;
}

template class MRGComponent<std::int64_t> ;
template class MRGComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif

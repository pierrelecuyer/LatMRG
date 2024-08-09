#ifndef	LATMRG_MRGCOMPONENT_H
#define	LATMRG_MRGCOMPONENT_H

#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactorization.h"
// #include "latmrg/Modulus.h"
#include <NTL/mat_poly_ZZ.h>
#include <string>
#include "Primitivity.h"

namespace LatMRG {

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

private:
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;

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
     * The names `filem1` and `filer` specify the file where the factors of
     * \f$m-1\f$ and \f$r\f$, are read or written when the type
     * <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt> is used.
     */
    MRGComponent(const Int &m, int k, DecompType decom1, const char *filem1,
            DecompType decor, const char *filer);

    /**
     * Same as the previous constructor with \f$m = b^e + c\f$.
     */
    MRGComponent(Int b, int e, Int c, int k, DecompType decom1,
            const char *filem1, DecompType decor, const char *filer);

    /**
     * Constructor similar to the above, except that the modulus of
     * congruence \f$m\f$ is inside the object `modul`.
     *               Should we keep this?      ************  ???
     */
//    MRGComponent(Modulus<Int> &modul, int k, DecompType decom1,
//            const char *filem1, DecompType decor, const char *filer);

    /**
     * Destructor.
     */
    ~MRGComponent();

    /**
     * Copy constructor.
     */
    MRGComponent(const MRGComponent<Int> &comp);

    /**
     * Assignment operator, assigns `comp` to this object.
     */
    MRGComponent& operator=(const MRGComponent<Int> &comp);

    /**
     * Sets the vector of multipliers to `aa`, with `aa[j]` containing \f$a_j\f$.
     * The order \f$k\f$ of the MRG is set equal to the length of this vector, plus 1.
     */
    void setaa(const IntVec &aa);

    /**
     * Returns the vector of multipliers of the recurrence, in same format as `aa` above.
     */
    IntVec getaa() const {
        return m_a;
    }

    /**
     * Returns the matrix multiplier for this MRG so it can be viewed as a MMRG.
     *                                                                Useful   ******  ??
     */
    IntMat getA() const {
        return m_A;     // Must be constructed!
    }

    /**
     * Sets the matrix of the recurrence to \f$A\f$.      NOT HERE! *******  ??
     * */
    // void setA(const IntMat &A);

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
    int getOrder() const {
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
    int& gete() {
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
     *            No underscore.  *****************
    void set_type(GenType type) {
        m_type = type;
    }
     * Gets the type of this component.
     *
    GenType get_type() {
        return m_type;
    }

    /**
     * The modulo of the MWC generator if we study one. This is because the
     * we check MWC period with module.m as the modulo of the equivalentn LCG.
     * */
    // Int m_MWCb;    // Does no belong here!     ************

    //    Does not belong here, I think.
    // void setPower2(std::vector<IntVec> &coeffs);

private:

    /**
     * This is called by the constructor, with the same arguments.
     */
    void init(const Int &m, int k, DecompType decom1, const char *filem1,
            DecompType decor, const char *filer);

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
    int m_m;

    /**
     * The order \f$k\f$ of the recurrence.
     */
    int m_k;

    /**
     * The multipliers \f$a_i\f$ of the recurrence, \f$i = 1, ..., k\f$.
     */
    IntVec m_a;

    /**
     * The generator matrix \f$A\f$ of the recurrence for MMRG   ******  ???
     */
    IntMat m_A;    // Want it here?

    /**
     * Basis `b` with `b^e+c = m`.
     * */
    Int m_b;

    /**
     * Exponent `e` with `b^e+c = m`.
     * */
    int m_e;

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


// Main constructor.
template<typename Int>
MRGComponent<Int>::MRGComponent(const Int &m, int k, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer) {
    m_b = m;
    m_e = 1;
    m_c = Int(0);
    init(m, k, decom1, filem1, decor, filer);
}

//===========================================================================

template<typename Int>
MRGComponent<Int>::MRGComponent(Int b, int e, Int c, int k, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer) {
    Int m = NTL::power(b, e) + c;
    m_b = b;
    m_e = e;
    m_c = c;
    init(m, k, decom1, filem1, decor, filer);
}

//===========================================================================

// Copy constructor.
template<typename Int>
MRGComponent<Int>::MRGComponent(const MRGComponent<Int> &lat) :
        ifm1(lat.ifm1), ifr(lat.ifr), m_k(lat.m_k) {
    m_m = lat.m_m;
    nj = lat.nj;
    rho = lat.rho;
    //   a.kill();
    m_a.resize(m_k);
    m_a = lat.m_a;
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

// Assignment.
template<typename Int>
MRGComponent<Int>& MRGComponent<Int>::operator=(const MRGComponent<Int> &lat) {
    if (this != &lat) {
        m_k = lat.m_k;
        m_m = lat.m_m;
        nj = lat.nj;
        rho = lat.rho;
        //    a.kill();
        m_a.resize(m_k);
        m_a = lat.m_a;
        //     orbitSeed.kill();
        orbitSeed.resize(m_k);
        orbitSeed = lat.orbitSeed;
        //      ifm1 = lat.ifm1;
        //      ifr = lat.ifr;
    }
    // m_type = lat.m_type;
    return *this;
}

//============================================================================


template<typename Int>
void MRGComponent<Int>::init(const Int &m0, int k0, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer) {
    FlexModInt<Int>::setModulus(m0);
    // module.init(m0);
    m_m = m0;
    m_k = k0;
    m_a.resize(m_k);
    orbitSeed.resize(m_k);

    Int m1;
    m1 = getModulus() - 1;
    ifm1.setNumber(m1);
    ifm2.setNumber(NTL::ZZ(m1));
    if (m_k == 1) {
        factor.setNumber(getModulus());
        factor.factorize();
    }
    if (decom1 != NO_DECOMP) {
        if (decom1 == DECOMP_READ)
            ifm1.read(filem1);
        else
            ifm1.factorize();
        if (decom1 == DECOMP_WRITE) {
            std::ofstream fout(filem1);
            fout << ifm1.toString();
        }
        {
            std::ofstream fout("dummy");
            fout << ifm1.toString();
        }
        ifm2.read("dummy");
        remove("dummy");
        ifm1.calcInvFactors();
        ifm2.calcInvFactors();
    }
    Int r;
    r = (NTL::power(m0, m_k) - 1) / (m0 - 1);
    ifr.setNumber(r);
    ifr2.setNumber(NTL::ZZ(r));

    if (decor != NO_DECOMP) {
        if (decor == DECOMP_READ)
            ifr.read(filer);
        else if (decor == DECOMP_PRIME)
            ifr.setStatus(LatticeTester::PRIME);
        else
            ifr.factorize();
        if (decor == DECOMP_WRITE) {
            std::ofstream fout(filer);
            fout << ifr.toString();
        }
        {
            std::ofstream fout("dummy");
            fout << ifr.toString();
        }
        ifr2.read("dummy");
        remove("dummy");
        ifr.calcInvFactors();
        ifr2.calcInvFactors();
    }

}

//===========================================================================

/**
template<typename Int>
MRGComponent<Int>::MRGComponent(Modulus<Int> &modu, int k, DecompType decom1,
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
    m_a.SetLength(aa.length());
    m_a = aa;
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

//  This one must be rewritten!       ****************
//  Maybe just replaced by a single call to `isPrimitive`, and that's all!

template<typename Int>
bool MRGComponent<Int>::maxPeriod(const IntVec &aa) {
    Primitivity<Int>::setM(getM());
    m_a = aa;      // The m_a is changed and this is not said in the doc!  ****
    Primitivity<Int>::reverse(m_a, m_k, 2);
    Primitivity<Int>::setF(m_a);
    Primitivity<Int> pol;
    PrimitiveInt<Int> privfm(ifm1, getM(), 1);
    return pol.isPrimitive(privfm, ifr);
}

//===========================================================================

//  This one is for an LCG, does not belong here.     ************
/*
template<typename Int>
bool MRGComponent<Int>::maxPeriod(const Int &a) {
    auto list = factor.getFactorList();
    for (auto iter = list.begin(); iter != list.end(); iter++) {
        if ((a0 - Int(1)) % (*iter).getFactor() != 0)
            return false;
    }
    if (getM() % 4 == 0 && (a0 - Int(1)) % 4 != 0)
        return false;

    return true;
}
*/

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
    for (int i = 0; i < m_k; i++) {
        for (int j = 0; j < m_k; j++) {
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
bool MRGComponent<Int>::maxPeriod23(const IntVec &a0) {
    FlexModInt<Int>::setModulus(getM());
    m_a = a0;
    PrimitivePoly<Int>::reverse(m_a, m_k, 2);
    PrimitivePoly<Int>::setF(m_a);
    PrimitivePoly<Int> pol;
    // La condition 1 a déjà été vérifiée.
    return pol.isPrimitive23(ifr);
}

//===========================================================================

template<typename Int>
std::string MRGComponent<Int>::toString() {
    std::ostringstream os;
    os << "MRGComponent:";
    Int mm = getM();
    os << "\n   m = " << mm << " = " << m_b << "^" << m_e << " + " << m_r;
    os << "\n   k = " << m_k;
    os << "\n   a = ";
    std::string str(os.str());
    std::string s2 = LatticeTester::toString(m_a, m_k);
    str += s2;
    str += "\n";
    return str;
}

extern template class MRGComponent<std::int64_t> ;
extern template class MRGComponent<NTL::ZZ> ;

} // End namespace LatMRG
#endif

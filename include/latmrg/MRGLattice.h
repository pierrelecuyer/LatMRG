#ifndef LATMRG_MRGLATTICE_H
#define LATMRG_MRGLATTICE_H

#include "latticetester/EnumTypes.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"

#include "latmrg/EnumTypes.h"
#include "latmrg/MRGComponent.h"

#include <string>

namespace LatMRG {

/**
 * This subclass of `IntLatticeExt` constructs the lattices associated with
 * a multiple recursive linear congruential generators (MRG), which follows
 * the recurrence
 * \f[
 *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mod m
 * \f]
 * where \f$m\f$ is the modulus, \f$k\f$  is the order, and  \f$a_1,\dots,a_k\f$
 * are the coefficients of the recurrence, also called the multipliers.
 * These coefficients are given in a vector `aa` where \f$a_j\f$ is in `aa[j-1]`.
 * The modulus \f$m\f$ and the order \f$k\f$ must be given to the constructor,
 * together with the maximal dimension `maxDim` for which the basis may be constructed.
 * We can then use the same `MRGLattice` objects for several different vectors
 * of coefficients by using the `setCoeff` and `buildBasis` or `buildDualBasis`
 * functions repeatedly. This is useful in particular when searching for good
 * vectors of multipliers for given values of \f$m\f$ and \f$k\f$.
 *
 * Only the `FULL` and `PRIMEPOWER` lattice types are currently implemented.
 * The others were partially implemented in older versions of the software.
 *
 * QUESTIONS / NOTES:
 */

template<typename Int, typename Real>
class MRGLattice: public LatticeTester::IntLatticeExt<Int, Real> {

public:

    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;

    /**
     * Constructs and `MRGLattice` object with modulus \f$m\f$, order \f$k\f$,
     * vector of multipliers \f$aa\f$, maximal dimension `MaxDim`,
     * lattice type `latType`, and using norm `norm`.
     * The coefficient \f$a_j\f$ must be given in `aa[j-1]`.
     * The `LatticeType` parameter is defined in `EnumTypes`, it can be:
     * FULL, RECURRENT, ORBIT, PRIMEPOWER (but only FULL and PRIMEPOWER is currently available).
     * The variable `withPrimal` indicates if the primal basis will be maintained or not,
     * and `withDual` indicates if the dual basis will be maintained or not.
     * This constructor does not build the basis, to leave
     * more flexibility in the dimension when doing so.
     * The internal vectors and matrices will be created with `maxDim` dimensions.

     * ****  indices in vectors and matrices vary from dimension 1 to `maxDim`. *****  CHANGE THIS?
     */
    MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, int64_t k = 1,
            bool withPrimal = false, bool withDual = false,
            LatticeType latType = FULL, LatticeTester::NormType norm =
                    LatticeTester::L2NORM);

    /**
     * Alternative constructor for an LCG. In this case, `k=1` and `a` is
     * a single integer.
     */
    MRGLattice(const Int &m, const Int &a, int64_t maxDim, bool withPrimal = false,
            bool withDual = false, LatticeType lattype = FULL,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

    /**
     * This constructor does not take a vector of coefficients.
     */
    MRGLattice(const Int &m, int64_t maxDim, int64_t k = 1, bool withPrimal = false,
            bool withDual = false, LatticeType lattype = FULL,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

    /**
     * Copy constructor.
     */
    MRGLattice(const MRGLattice<Int, Real> &lat);

    /**
     * Assigns `lat` to this object (deep copy).
     */
    MRGLattice<Int, Real>& operator=(const MRGLattice<Int, Real> &Lat);

    /**
     * Destructor.
     */
    ~MRGLattice();

    /**
     * Cleans and releases the memory used by this object.
     */
    void kill();

    /**
     * Sets the vector of multipliers to `aa`.  This vector must have `k` dimensions.
     */
    void setaa(const IntVec &aa);

    /**
     * Assumes that `k=1` and sets this object to an LCG with multiplier `a`.
     */
    void seta(const Int &a);

    /**
     * Sets `m_power2` to true and sets `m_pow2_exp to `coeffs`.  ******  ????
     */
    void setPower2(std::vector<IntVec> &coeffs);

    /**
     * Builds a basis in `dim` dimensions. This `dim` must not exceed `maxDim`.
     * This initial primal basis will be upper triangular and built as explained
     * in the LatMRG guide.
     * This function can be called only when `withPrimal` is set to true.
     * If `withDual` is true, it also builds an m-dual lower-triangular basis.
     */
    virtual void buildBasis(int64_t dim);

    /**
     * Builds only an m-dual lower triangular basis (and not the primal) directly
     * in `dim` dimensions, as explained in the LatMRG guide.
     * This `dim` must not exceed `maxDim`.
     */
    void buildDualBasis(int64_t dim);

    /**
     * Increases the current dimension of the primal basis by 1 and updates the basis.
     * This function can be called only when `withPrimal` is set to true.
     * If `withDual`, it also increases the m-dual basis and makes it the m-dual of the primal basis.
     * The new increased dimension must not exceed `maxDim`.
     */
    void incDimBasis();

    /**
     * Increases the current dimension of only the m-dual basis by 1.
     * The primal basis is left unchanged (not updated).
     * The new increased dimension must not exceed `maxDim`.
     * This function uses the direct method given in the LatMRG guide.
     */
    void incDimDualBasis();

    /**
     * This function computes a basis for the projection of this lattice over the
     * coordinates in `proj`, and returns it in the `projLattice` object.
     * The implementation exploits the lattice structure of the MRG as explained in
     * the guide. If `withDual` is true, it also computes an m-dual basis for the
     * projection and stores it in `projLattice`.
     * In general, a set of generating vectors for the primal lattice is obtained
     * by picking the columns whose indices are in `proj`, and reducing the resulting
     * matrix to a basis using the tools from `BasisConstruction`.
     * If we only want a primal basis, we should use LLL for the reduction.
     * If we also want the m-dual, then it is better to compute an upper-triangular
     * basis for the primal and the corresponding lower-triangular m-dual basis.
     * This function does not require that a basis for the whole lattice was constructed before.
     */
    void buildProjection(IntLattice<Int, Real> *projLattice,
            const Coordinates &proj, double delta = 0.99) override;

    /**
     * Returns a reference to the vector of multipliers of this MRG.
     */
    const IntVec& getCoeff() const {
        return m_aCoeff;
    }

    /**
     * Returns the vector of multipliers \f$\ba \f$ as a string.
     */
    virtual std::string toStringCoeff() const;

    /**
     * Returns the vector of multipliers \f$A\f$ as a string.   ????????
     */
    std::string toString() const override;

protected:

    /**
     * Initializes a square matrix of order \f$k\f$. This initial matrix
     * contains a system of generators for the given group of states.
     */
    void initStates();

    /**
     * Initializes some of the local variables.
     */
    void init();

    /**
     * Initializes this object when the lattice type is `ORBIT`.
     */
    // void initOrbit();
    // void insertion(IntVec &Sta);
    // void lemme2(IntVec &Sta);

    /**
     * For debugging purposes.
     */
    void trace(char *msg, int64_t d);

    /**
     * Increments the basis by 1 in case of non-lacunary indices.
     */
    virtual void incDimBasis();

    /**
     * The coefficients of the recurrence.
     */
    IntVec m_aCoeff;

    /**
     * Indicates which lattice or sublattice is analyzed.
     */
    LatticeType m_latType;

    /**
     * Matrix that contains the vectors that can be used to generate the
     * basis for an arbitrary dimension. This matrix is of order k and if
     * we want to build the full lattice, this matrix is the identity
     * matrix. **Marc-Antoine** This matrix is different in some way that I
     * don't quite understand if we use lacunary indices.
     */
    IntMat m_sta;

    /**
     * When the flag `m_ip[i]` is `true`, the \f$i\f$-th diagonal
     * element of matrix `m_sta` is non-zero (modulo \f$m\f$) and
     * divides \f$m\f$. Otherwise (when `m_ip[i]` is
     * `false`), the \f$i\f$-th line of matrix `m_sta` is
     * identically 0.
     */
    bool *m_ip;

private:

    /**
     * Indicates if this generator is built using coefficients that are
     * a sum or difference of a small number of powers of 2.
     */
    bool m_power2;

    /**
     * The powers of 2 used if this generator uses power of 2 coefficients.
     */
    std::vector<IntVec> m_pow2_exp;

};
// End class declaration


//===========================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &a, int64_t maxDim,
        int64_t k, LatticeType lat, LatticeTester::NormType norm) :
        LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(m, k, maxDim,
                true, norm) {
    m_latType = lat;
    m_lacunaryFlag = false;
    m_ip = 0;
    init();

    for (int64_t i = 0; i < this->m_order; i++)
        m_aCoeff[i] = a[i + 1];
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const Int &a, int64_t maxDim,
        LatticeType lat, LatticeTester::NormType norm) :
        LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(m, 1, maxDim,
                true, norm) {
    m_latType = lat;
    m_lacunaryFlag = false;
    m_ip = 0;
    init();
    m_aCoeff[0] = a;
}

//===========================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &a, int64_t maxDim,
        int64_t k, IntVec &I, LatticeType lat, LatticeTester::NormType norm) :
        LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(m, k, maxDim,
                true, norm), m_lac(I, maxDim), m_ip(0) {
    std::cout << 4 << "\n";
    m_latType = lat;
    m_lacunaryFlag = true;
    init();

    for (int64_t i = 0; i < this->m_order; i++)
        m_aCoeff[i] = a[i + 1];
}

//===========================================================================

// Copy constructor, makes a deep copy.
template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
        LatticeTester::IntLatticeExt<Int, Real>::IntLatticeExt(lat.m_modulo,
                lat.m_order, lat.getDim(), lat.withDual(), lat.getNorm()), m_lac(
                lat.m_lac) {
    m_latType = lat.m_latType;
    m_lacunaryFlag = lat.m_lacunaryFlag;
    init();

    for (int64_t i = 0; i < this->m_order; i++)
        m_aCoeff[i] = lat.m_aCoeff[i];

    m_power2 = lat.m_power2;
    if (m_power2) {
        m_pow2_exp.resize(lat.m_pow2_exp.size());
        for (unsigned int64_t i = 0; i < lat.m_pow2_exp.size(); i++) {
            m_pow2_exp[i].resize(lat.m_pow2_exp[i].length());
            for (int64_t j = 0; j < lat.m_pow2_exp[i].length(); j++)
                m_pow2_exp[i][j] = lat.m_pow2_exp[i][j];
        }
    }
}

//===========================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(
        const MRGLattice<Int, Real> &lat) {
    if (this == &lat)
        return *this;
    kill();
    m_latType = lat.m_latType;
    m_lacunaryFlag = lat.m_lacunaryFlag;

    for (int64_t i = 0; i < this->m_order; i++)
        m_aCoeff[i] = lat.m_aCoeff[i];
    this->m_dim = lat.m_dim;
    this->copyBasis(lat);
    this->m_order = lat.m_order;
    this->m_ip = lat.m_ip;

    m_power2 = lat.m_power2;
    if (m_power2) {
        m_pow2_exp.resize(lat.m_pow2_exp.size());
        for (unsigned int64_t i = 0; i < lat.m_pow2_exp.size(); i++) {
            m_pow2_exp[i].resize(lat.m_pow2_exp[i].length());
            for (int64_t j = 0; j < lat.m_pow2_exp[i].length(); j++)
                m_pow2_exp[i][j] = lat.m_pow2_exp[i][j];
        }
    }
    //m_shift = lat.m_shift;
    return *this;
    //MyExit (1, " MRGLattice::operator= n'est pas terminé   " );
    //copy (lat);
    //return *this;
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::init() {
    m_xi.SetLength(this->m_order);
    m_aCoeff.SetLength(this->m_order);
    if (this->m_order > ORDERMAX) {
        m_ip = new bool[1];
        m_sta.SetDims(1, 1);
    } else {
        m_ip = new bool[this->m_order];
        m_sta.SetDims(this->m_order, this->m_order);
    }
    int64_t rmax = std::max(this->m_order, this->getDim());
    this->m_wSI.SetDims(rmax, this->getDim());

    m_power2 = false;
    if (m_latType == ORBIT)
        initOrbit();
}

//===========================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::~MRGLattice() {
    kill();
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::kill() {
    if (0 != m_ip)
        delete[] m_ip;
    m_ip = 0;
    m_xi.kill();
    m_aCoeff.kill();
    m_sta.kill();
}

//===========================================================================

template<typename Int, typename Real>
Int& MRGLattice<Int, Real>::getLac(int64_t j) {
    if (isLacunary() && j <= m_lac.getSize() && j > 0)
        return m_lac.getLac(j);
    throw std::out_of_range("MRGLattice::getLac");
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::setLac(const LatticeTester::Lacunary<Int> &lac) {
    m_lac = lac;
    m_lacunaryFlag = true;
}

//===========================================================================

template<typename Int, typename Real>
std::string MRGLattice<Int, Real>::toStringCoeff() const {
    std::ostringstream out;
    //out << "[ ";
    for (int64_t i = 0; i < this->m_order; i++)
        out << m_aCoeff[i] << "  ";
    //out << "]";
    return out.str();
}

//===========================================================================

template<typename Int, typename Real>
std::string MRGLattice<Int, Real>::toString() const {
    std::ostringstream out;
    for (int64_t i = 0; i < this->m_order; i++) {
        out << "a_" << i + 1 << " = " << m_aCoeff[i];
        if (m_power2 && m_aCoeff[i] != 0) {
            out << " = ";
            Int tmp = m_aCoeff[i];
            for (int64_t j = 0; j < m_pow2_exp[i].length(); j++) {
                Int tmp2, tmp3 = m_pow2_exp[i][j];
                NTL::power2(tmp2, tmp3);
                if (tmp > 0) {
                    if (j != 0)
                        out << " + ";
                    out << "2^" << m_pow2_exp[i][j];
                    tmp -= tmp2;
                } else if (tmp < 0) {
                    if (j != 0)
                        out << " - ";
                    else
                        out << "-";
                    out << "2^" << m_pow2_exp[i][j];
                    tmp += tmp2;
                } else
                    out << 0;
            }
        }
        out << "\n";
    }
    return out.str();
}


//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   // La base est construite en dimension d.

    initStates();

    int64_t dk = d;
    if (dk > this->m_order)
        dk = this->m_order;
    this->m_basis.SetDims(dk, dk);
    this->m_dualbasis.SetDims(dk, dk);
    this->m_vecNorm.SetLength(dk);
    this->m_dualvecNorm.SetLength(dk);
    this->setDim(dk);
    this->setNegativeNorm();
    this->setDualNegativeNorm();

    int64_t i, j;
    for (i = 0; i < dk; i++) {
        for (j = 0; j < dk; j++) {
            if (i == j) {
                this->m_basis[i][j] = Int(1);
                this->m_dualbasis[i][j] = this->m_modulo;
            } else {
                this->m_basis[i][j] = Int(0);
                this->m_dualbasis[i][j] = Int(0);
            }
        }
    }

    if (d > this->m_order) {
        for (i = this->m_order; i < d; i++)
            incDimBasis();
    }
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDim() {
    if (m_lacunaryFlag) {
        incDimLaBasis(this->getDim());
    } else {
        incDimBasis();
    }
    //   write (1);
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
    LatticeTester::IntLatticeExt<Int, Real>::incDim();

    const int64_t dim = this->getDim();
    this->m_vSI.resize(dim, dim);
    this->m_wSI.resize(dim, dim);
    this->m_vecNorm.resize(dim);
    this->m_dualvecNorm.resize(dim);

    for (int64_t i = 0; i < dim; i++) {
        NTL::clear(this->m_vSI[0][i]);
        for (int64_t j = 0; j < this->m_order; j++) {
            NTL::conv(this->m_t1, this->m_basis[i][dim - j - 2]);
            this->m_t1 = this->m_t1 * m_aCoeff[j];
            this->m_vSI[0][i] = this->m_vSI[0][i] + this->m_t1;
        }
        LatticeTester::Modulo(this->m_vSI[0][i], this->m_modulo,
                this->m_vSI[0][i]);
        this->m_basis[i][dim - 1] = this->m_vSI[0][i];
    }

    for (int64_t i = 0; i < dim; i++)
        this->m_basis[dim - 1][i] = 0;
    this->m_basis[dim - 1][dim - 1] = this->m_modulo;

    for (int64_t i = 0; i < dim - 1; i++) {
        this->m_dualbasis[i][dim - 1] = 0;
        this->m_dualbasis[dim - 1][i] = -this->m_vSI[0][i];
    }
    this->m_dualbasis[dim - 1][dim - 1] = 1;

    /*
     * This old version used the dual to compute the next vector in the dual
     * even when it was not needed. This caused problems because we could not
     * increase the lattice dimension after performing computations on the
     * dual because this would not add the correct column to the dual.
     * */
    // for (int64_t j = 0; j < dim-1; j++) {
    //   NTL::clear (this->m_t1);
    //   for (int64_t i = 0; i < dim-1; i++) {
    //     this->m_t2 = this->m_dualbasis[i][j];
    //     this->m_t2 *= this->m_vSI[0][i];
    //     this->m_t1 -= this->m_t2;
    //   }
    //   LatticeTester::Quotient (this->m_t1, this->m_modulo, this->m_t1);
    //   this->m_dualbasis[dim-1][j] = this->m_t1;
    // }
    this->setNegativeNorm();
    this->setDualNegativeNorm();
    /*
     if (!checkDuality ())
     MyExit (1, "BUG");
     */
    // trace("=================================APRES incDimBasis", -10);
}

//===========================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(
        LatticeTester::IntLatticeExt<Int, Real> *lattice,
        const LatticeTester::Coordinates &proj) {
    // Ripping the vectors
    IntMat temp(this->m_dim, proj.size());
    IntMat temp2(proj.size(), proj.size());
    int64_t i = 0;
    for (auto iter = proj.begin(); iter != proj.end(); ++iter) {
        for (int64_t j = 0; j < this->m_dim; j++) {
            temp[j][i] = this->m_basis[j][*iter];
        }
        i++;
    }

    // Construct the basis for the projection, via LLL.
    LatticeTester::BasisConstruction<Int> constr;
    constr.LLLConstruction(temp);
    constr.DualConstruction(temp, temp2, this->m_modulo);
    for (int64_t i = 0; i < lattice->getDim(); i++) {
        for (int64_t j = 0; j < lattice->getDim(); j++) {
            lattice->getBasis()[i][j] = temp[i][j];
            lattice->getDualBasis()[i][j] = temp2[i][j];
        }
    }
}


//===========================================================================

// The functions initStates below works only for the cases FULL and PRIMEPOWER.
// The functions insertion, lemme2, initOrbit can be ignored for now.     ******
// They are used only for the ORBIT case, which is currently not working. ******

// Initialise la matrice carr�e Sta, dont l'ordre est �gal a�l'ordre k du générateur.
// Elle contient un système de générateurs pour le groupe d'états considérés.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::initStates() {
    IntVec statmp;
    statmp.resize(this->m_order);    // To store some variables
    int64_t maxDim = this->getDim();
    //clear (m_t2);

    if (m_latType == FULL || m_latType == PRIMEPOWER
            || (this->m_t1 == this->m_t2)) {
        // We want the full lattice that contains all the orbits.
        // m_sta is set to identity matrix
        for (int64_t i = 0; i < this->m_order; i++) {
            for (int64_t j = 0; j < this->m_order; j++) {
                if (i != j)
                    NTL::clear(m_sta[i][j]);
                else
                    NTL::set(m_sta[i][j]);
            }
            m_ip[i] = true;
        }
        double temp;
        NTL::conv(temp, this->m_modulo);
        double lgm2 = 2.0 * LatticeTester::Lg(temp);
        this->calcLgVolDual2(lgm2);

        for (int64_t r = this->m_order + 1; r <= maxDim; r++)
            this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
    }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::setPower2(std::vector<IntVec> &coeffs) {
    this->m_power2 = true;
    this->m_pow2_exp.resize(this->m_order);
    for (int64_t i = 0; i < this->m_order; i++) {
        this->m_pow2_exp[i].resize(coeffs[i].length());
        // Sorting the exponents
        for (int64_t j = 0; j < coeffs[i].length(); j++) {
            for (int64_t k = j; k < coeffs[i].length(); k++) {
                if (m_pow2_exp[i][j] < coeffs[i][k]) {
                    m_pow2_exp[i][j] = coeffs[i][k];
                    coeffs[i][k] = coeffs[i][j];
                    coeffs[i][j] = m_pow2_exp[i][j];
                }
            }
        }
    }
}

template class MRGLattice<std::int64_t, double> ;
template class MRGLattice<NTL::ZZ, double> ;
template class MRGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatMRG
#endif

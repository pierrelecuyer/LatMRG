#ifndef	LATMRG_MRGCOMPONENT_H
#define	LATMRG_MRGCOMPONENT_H

#include "latmrg/Const.h"
#include "latmrg/PolyPE.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Modulus.h"
#include "latmrg/IntPrimitivity.h"

#include <NTL/mat_poly_ZZ.h>

#include <string>

namespace LatMRG {

  /**
   * This class is used to to test for full period length of generators and also
   * to store a bit of information about generators without creating a lattice.
   * Each MRG component is defined by a
   * modulus \f$m\f$, an order \f$k\f$ and a list of files and factorizations
   * needed check the full period.
   *
   * This class can also store
   * - a vector of multipliers \f$a\f$,
   * where \f$a[i-1]\f$ represents \f$a_i\f$. This MRG satisfies the recurrence
   * \f[
   *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mbox{ mod } m
   * \f]
   * - or a matrix A in case of MMRG
   *
   * If it is used to check the period of another generator, the multipliers
   * stored might change.
   */
  template<typename Int>
    class MRGComponent {
      private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;
      public:

        /**
         * Constructor with modulus \f$m\f$ and order \f$k\f$. Arguments
         * `decom1` and `decor` refer to the prime factor decomposition of
         * \f$m-1\f$ and \f$r=(m^k-1)/(m-1)\f$, respectively. If `decor` equals
         * `DECOMP`, the constructor will factorize \f$r\f$. If `decor` equals
         * <tt>DECOMP_WRITE</tt>, the constructor will factorize \f$r\f$ and
         * write the prime factors to file `filer`. If `decor` equals
         * <tt>DECOMP_READ</tt>, the constructor will read the factors of
         * \f$r\f$ from file `filer`. If `decor` equals <tt>DECOMP_PRIME</tt>,
         * \f$r\f$ is assumed to be prime. Similar considerations apply to
         * `decom1` and `filem1` with respect to \f$m-1\f$.
         */
        MRGComponent (const Int & m, int k, DecompType decom1,
            const char *filem1, DecompType decor,
            const char *filer);

        /**
         * Same as the other constructors with `m=b^e+r`.
         * */
        MRGComponent (Int b, int e, Int r, int k, DecompType decom1,
            const char *filem1, DecompType decor,
            const char *filer);

        /**
         * Constructor similar to the above, except that the modulus of
         * congruence \f$m\f$ is inside the object `modul`.
         */
        MRGComponent (Modulus<Int> & modul, int k, DecompType decom1,
            const char *filem1, DecompType decor,
            const char *filer);

        /**
         * Destructor.
         */
        ~MRGComponent();

        /**
         * Copy constructor;
         */
        MRGComponent (const MRGComponent<Int> & comp);

        /**
         * Assignment operator.
         */
        MRGComponent<Int> & operator= (const MRGComponent<Int> & comp);

        /**
         * Sets the multipliers of the recurrence to \f$A\f$.
         */
        void setA (const IntVec & A);

        /**
         * Gets the multipliers of the recurrence to \f$A\f$.
         */
        IntVec getA() const {return m_a;}

        /**
         * Gets the matrix multipliers of the recurrence.
         * */
        IntMat getMatrix() const {return m_A;}

        /**
         * Sets the matrix of the recurrence to \f$A\f$.
         * */
        void setA(const IntMat& A);

        /**
         * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
         * period; returns `false` otherwise. Note that when calling this class
         * it is necessary that `a[i] = a_i`.
         */
        bool maxPeriod (const IntVec & A);

        /**
         * Returns `true` if \f$a\f$ makes for a full period LCG with non-zero
         * carry; returns `false` otherwise.
         *
         * This checks the 2 following conditions :
         * - If \f$q\f$ is prime and divides \f$m\f$, it must divide \f$a-1\f$.
         * - If 4 divides \f$m\f$, it must divide \f$a-1\f$.
         *
         * The user must choose the carry himself. Any carry relatively prime to
         * \f$m\f$ will give full period.
         */
        bool maxPeriod (const Int & a);

        /**
         * Returns `true` if coefficients in \f$A\f$ give a MMRG with maximal
         * period; returns `false` otherwise.
         */
        bool maxPeriod (const IntMat & A);

        /**
         * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
         * period; returns `false` otherwise. This method supposes that
         * condition 1 is `true` and tests only conditions 2 and 3. See method
         * `isPrimitive` of class `PolyPE` on page (FIXME: page#) of this
         * guide.
         */
        bool maxPeriod23 (const IntVec & A);

        /**
         * Returns the value of the modulus \f$m\f$ of the recurrence.
         */
        const Int getM() const { return module.m; }

        /**
         * Returns the value of the modulus \f$m\f$ of the recurrence.
         */
        int getK() const { return m_k; }

        /**
         * Returns the value of `b` with `b^e+r = m`.
         */
        Int getB() const { return m_b; }

        /**
         * Returns the value of `e` with `b^e+r = m`.
         */
        int getE() const { return m_e; }

        /**
         * Returns the value of `r` with `b^e+r = m`.
         */
        Int getR() const { return m_r; }

        /**
         * Returns const reference of orbitSeed
         * */
        IntVec& getOrbitSeed() { return orbitSeed;}

        /**
         * Returns this object as a string.
         */
        std::string toString ();

        /**
         * Sets the type of this component.
         * */
        void set_type(GenType type) {m_type = type;}

        /**
         * Gets the type of this component.
         * */
        GenType get_type() {return m_type;}

      private:

        /**
         * Does the same as the constructor above with similar arguments.
         */
        void init (const Int & m, int k, DecompType decom1,
            const char *filem1, DecompType decor,
            const char *filer);

        /**
         * The type of generator this stores. Should be MRG, MMRG or MWC. The
         * default value is LCG. When the value is LCG, this means the object is
         * not used to represent a generator, but to compute period length.
         * */
        GenType m_type = LCG;

        /**
         * The prime factor decomposition of \f$m-1\f$.
         */
        IntFactorization<Int> ifm1;

        /**
         * The prime factor decomposition of \f$m\f$.
         * Used to compute the full period of a LCG with a carry.
         */
        IntFactorization<Int> factor;

        /**
         * The prime factor decomposition of \f$m-1\f$.
         * This is used when computing the full period of a matrix generator.
         * This is because we use NTL to compute the characteristic polynomial
         * of the matrix since it is already implemented.
         */
        IntFactorization<NTL::ZZ> ifm2;

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
        IntFactorization<NTL::ZZ> ifr2;

        /**
         * The modulus \f$m\f$ of the recurrence.
         */
        Modulus<Int> module;

        /**
         * The order \f$k\f$ of the recurrence.
         */
        int m_k;

        /**
         * The multipliers \f$a_i\f$ of the recurrence, \f$i = 1, …, k\f$.
         */
        IntVec m_a;

        /**
         * The generator matrix \f$A\f$ of the recurrence for MMRG.
         */
        IntMat m_A;

        /**
         * Basis `b` with `b^e+r = m`.
         * */
        Int m_b;

        /**
         * Exponent `e` with `b^e+r = m`.
         * */
        int m_e;

        /**
         * Rest `r` with `b^e+r = m`.
         * */
        Int m_r;

        /**
         * The length of the period \f$\rho\f$ for this MRG. For now, the
         * recurrence is assumed to correspond to a primitive polynomial and
         * `rho` is calculated as
         * \f[
         * \rho= m^k - 1
         * \f]
         * This value is calculated by `MRGLatticeFactory` and stored here for
         * simplicity.
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
        Int nj;

        /**
         * Contains the starting state of the component for the case when the
         * lattice type is `ORBIT`. It is made of \f$k\f$ numbers.
         */
        IntVec orbitSeed;

    }; // End class declaration

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const MRGComponent<Int> & lat) :
      ifm1(lat.ifm1), ifr(lat.ifr), m_k(lat.m_k)
  {
    module = lat.module;
    nj = lat.nj;
    rho = lat.rho;
    //   a.kill();
    m_a.resize(m_k);
    LatticeTester::CopyVect(m_a, lat.m_a, m_k);
    //   orbitSeed.kill();
    orbitSeed.resize(m_k);
    LatticeTester::CopyVect(orbitSeed, lat.orbitSeed, m_k);

    m_type = lat.m_type;
    m_b = lat.m_b;
    m_e = lat.m_e;
    m_r = lat.m_r;
  }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int> & MRGComponent<Int>::operator=
    (const MRGComponent<Int> & lat)
    {
      if (this != &lat) {
        m_k = lat.m_k;
        module = lat.module;
        nj = lat.nj;
        rho = lat.rho;
        //    a.kill();
        m_a.resize(m_k);
        LatticeTester::CopyVect(m_a, lat.m_a, m_k);
        //     orbitSeed.kill();
        orbitSeed.resize(m_k);
        LatticeTester::CopyVect(orbitSeed, lat.orbitSeed, m_k);
        //      ifm1 = lat.ifm1;
        //      ifr = lat.ifr;
      }
      m_type = lat.m_type;
      return *this;
    }

  //============================================================================

  template<typename Int>
    void MRGComponent<Int>::init (const Int & m0, int k0, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer)
    {
      PolyPE<Int>::setM(m0);
      module.init(m0);
      m_k = k0;
      m_a.resize(m_k);
      orbitSeed.resize(m_k);

      Int m1;
      m1 = getM() - 1;
      ifm1.setNumber (m1);
      ifm2.setNumber (NTL::ZZ(m1));

      if (m_k == 1) {
        factor.setNumber(getM());
        factor.factorize();
      }

      if (decom1 != NO_DECOMP) {
        if (decom1 == DECOMP_READ) 
          ifm1.read (filem1);
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
          ifr.read (filer);
        else if (decor == DECOMP_PRIME) 
          ifr.setStatus (LatticeTester::PRIME);
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

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const Int & m, int k, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer)
    {
      m_b = m;
      m_e = 1;
      m_r = Int(0);
      init (m, k, decom1, filem1, decor, filer);
    }

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (Int b, int e, Int r, int k, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer)
    {
      Int m = NTL::power(b,e) + r;
      m_b = b;
      m_e = e;
      m_r = r;
      init (m, k, decom1, filem1, decor, filer);
    }

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (Modulus<Int> & modu, int k,
        DecompType decom1, const char *filem1, DecompType decor,
        const char *filer)
    {
      init (modu.m, k, decom1, filem1, decor, filer);

      LatticeTester::PrimeType status = LatticeTester::IntFactor<Int>::isPrime (
          modu.m, 100);
      if (status == LatticeTester::PRIME
          || LatticeTester::PROB_PRIME == status) {
        modu.primeF = true;
      } else {
        std::cout << " WARNING:  m is NOT prime" << std::endl;
        modu.primeF = false;
      }
      module = modu;
    }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::~MRGComponent()
    {
      //a.kill();
      //orbitSeed.kill();
    }


  //===========================================================================

  template<typename Int>
    void MRGComponent<Int>::setA (const IntVec & b)
    {
      m_a.SetLength(b.length());
      m_a = b;
      //  LatticeTester::CopyVect(b, a, k);
    }

  //===========================================================================

  template<typename Int>
    void MRGComponent<Int>::setA (const IntMat & b)
    {
      m_A.SetDims(b.NumRows(), b.NumCols());
      m_A = b;
      //  LatticeTester::CopyVect(b, a, k);
    }

  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod (const IntVec & a0)
    {
      PolyPE<Int>::setM(getM());
      m_a = a0;
      PolyPE<Int>::reverse (m_a, m_k, 2);
      PolyPE<Int>::setF(m_a);
      PolyPE<Int> pol;
      IntPrimitivity<Int> privfm(ifm1, getM(), 1);
      return pol.isPrimitive(privfm, ifr);
    }

  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod (const Int& a0)
    {
      auto list = factor.getFactorList();
      for (auto iter = list.begin(); iter != list.end(); iter++) {
        if ((a0-Int(1))%(*iter).getFactor() != 0) return false;
      }
      if (getM() % 4 == 0 && (a0-Int(1))%4 != 0) return false;

      return true;
    }

  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod (const IntMat & a0)
    {
      // Conerting everything to NTL::ZZ to ease the characteristic polynomial
      // computation
      PolyPE<NTL::ZZ>::setM(NTL::ZZ(getM()));
      NTL::ZZX poly;
      NTL::matrix<NTL::ZZ> mat(m_k, m_k);
      for (int i = 0; i<m_k; i++) {
        for (int j = 0; j<m_k; j++) {
          mat[i][j] = NTL::ZZ(a0[i][j]);
        }
      }
      // Characteristic polynomial computation
      NTL::CharPoly(poly, mat);
      // Copying the polynomial to a vector
      NTL::vector<NTL::ZZ> vec(NTL::VectorCopy(poly, m_k+1));
      PolyPE<NTL::ZZ>::setF(vec);
      PolyPE<NTL::ZZ> pol;
      IntPrimitivity<NTL::ZZ> privfm(ifm2, NTL::ZZ(getM()), 1);
      return pol.isPrimitive(privfm, ifr2);
    }


  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod23 (const IntVec & a0)
    {
      PolyPE<Int>::setM(getM());
      m_a = a0;
      PolyPE<Int>::reverse (m_a, m_k, 2);
      PolyPE<Int>::setF(m_a);
      PolyPE<Int> pol;
      // La condition 1 a déjà été vérifiée dans SeekMain
      return pol.isPrimitive(ifr);
    }


  //===========================================================================

  template<typename Int>
    std::string MRGComponent<Int>::toString ()
    {
      std::ostringstream os;
      os << "MRGComponent:";
      Int mm = getM();
      os << "\n   m = " << mm << " = " << m_b << "^" << m_e << " + " << m_r;
      os << "\n   k = " << m_k;
      os << "\n   a = ";
      std::string str (os.str ());
      std::string s2 = LatticeTester::toString(m_a, m_k);
      str += s2;
      str += "\n";
      return str;
    }

  extern template class MRGComponent<std::int64_t>;
  extern template class MRGComponent<NTL::ZZ>;

} // End namespace LatMRG
#endif

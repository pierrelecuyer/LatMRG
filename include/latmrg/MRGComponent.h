#ifndef	MRGCOMPONENT_H
#define	MRGCOMPONENT_H
#include "latticetester/Types.h"
#include "latticetester/Const.h"

#include "latmrg/Const.h"
#include "latmrg/PolyPE.h"
#include "latmrg/IntFactorization.h"
#include "latmrg/Modulus.h"

#include <string>


namespace LatMRG {

  /**
   * This class is used to implement a MRG component in a combined MRG. It
   * exists in order to avoid creating numerous relatively heavy `MRGLattice`
   * objects to represent MRG components. Each MRG component is defined by a
   * modulus \f$m\f$, an order \f$k\f$ and 
   * - either a vector of multipliers \f$a\f$,
   * where \f$a[i]\f$ represents \f$a_i\f$. This MRG satisfies the recurrence
   * \f[
   *   x_n = (a_1 x_{n-1} + \cdots+ a_k x_{n-k}) \mbox{ mod } m
   * \f]
   * - or a matrix A in case of MMRG
   */
  template<typename Int>
    class MRGComponent {
      public:

        /**
         * Constructor with modulus \f$m\f$, vector \f$a\f$ and order \f$k\f$.
         */
        MRGComponent (const Int & m, const MVect & a, int k);


        /**
         * Constructor for MMRG with modulus \f$m\f$, Matrix \f$A\f$ and
         * order \f$k\f$.
         */
        MRGComponent (const Int & m, const MMat & A, int k);

        /**
         * Constructor for MMRG with modulus \f$m\f$, Matrix \f$A\f$ and
         * order \f$k\f$.
         */
        MRGComponent (long b, long e, long c, const MMat & A, int k);


        /**
         * Constructor with modulus \f$m=b^e + c\f$, vector \f$a\f$ and order
         * \f$k\f$.
         */
        MRGComponent (long b, long e, long c, const MVect & a, int k);

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
        void setA (const MVect & A);

        /**
         * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
         * period; returns `false` otherwise.
         */
        bool maxPeriod (const MVect & A);

        /**
         * Returns `true` if coefficients \f$A\f$ give a MRG with maximal
         * period; returns `false` otherwise. This method supposes that
         * condition 1 is `true` and tests only conditions 2 and 3. See method
         * `isPrimitive` of class `PolyPE` on page (FIXME: page#) of this
         * guide.
         */
        bool maxPeriod23 (const MVect & A);

        /**
         * The prime factor decomposition of \f$m-1\f$.
         */
        IntFactorization<Int> ifm1;

        /**
         * The prime factor decomposition of \f$r=(m^k-1)/(m-1)\f$, where
         * \f$k\f$ is the order of the recurrence.
         */
        IntFactorization<Int> ifr;

        /**
         * The modulus \f$m\f$ of the recurrence.
         */
        Modulus<Int> module;

        /**
         * Returns the value of the modulus \f$m\f$ of the recurrence.
         */
        Int getM() { return module.m; }

        /**
         * The order \f$k\f$ of the recurrence.
         */
        int k;

        /**
         * The multipliers \f$a_i\f$ of the recurrence, \f$i = 1, …, k\f$.
         */
        MVect a;

        /**
         * The generator matrix \f$A\f$ of the recurrence for MMRG.
         */
        MMat A;

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
        MVect orbitSeed;

        /**
         * Returns this object as a string.
         */
        std::string toString ();
      private:

        /**
         * Does the same as the constructor above with similar arguments.
         */
        void init (const Int & m, int k, DecompType decom1,
            const char *filem1, DecompType decor,
            const char *filer);
    }; // End class declaration

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const Int & m, const MVect & a0, int k0)
    {
      PolyPE::setM(m);
      module.init(m);
      k = k0;
      a.resize(k);
      LatticeTester::CopyVect(a, a0, k);
      orbitSeed.resize(k);
    }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const Int & m, const MMat & A0, int k0)
    {
      PolyPE::setM(m);
      module.init(m);
      k = k0;
      A.resize(k, k);
      LatticeTester::CopyMatr(A, A0, k);
      orbitSeed.resize(k);
    }

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (long p, long e, long c, const MMat & A0, 
        int k0)
    {
      module.init(p, e, c);
      PolyPE::setM(getM());
      k = k0;
      A.resize(k, k);
      LatticeTester::CopyMatr(A, A0, k);
      orbitSeed.resize(k);
    }

  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (long p, long e, long c, const MVect & a0, 
        int k0)
    {
      module.init(p, e, c);
      PolyPE::setM(getM());
      k = k0;
      a.resize(k);
      LatticeTester::CopyVect(a, a0, k);
      orbitSeed.resize(k);
    }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const MRGComponent<Int> & lat) : 
      k(lat.k)//, ifm1(lat.ifm1), ifr(lat.ifr)
  {
    module = lat.module;
    nj = lat.nj;
    rho = lat.rho;
    //   a.kill();
    a.resize(k + 1);
    LatticeTester::CopyVect(a, lat.a, k);
    //   orbitSeed.kill();
    orbitSeed.resize(k + 1);
    LatticeTester::CopyVect(orbitSeed, lat.orbitSeed, k);
  }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int> & MRGComponent<Int>::operator= 
    (const MRGComponent<Int> & lat)
    {
      if (this != &lat) {
        k = lat.k;
        module = lat.module;
        nj = lat.nj;
        rho = lat.rho;
        //    a.kill();
        a.resize(k + 1);
        LatticeTester::CopyVect(a, lat.a, k);
        //     orbitSeed.kill();
        orbitSeed.resize(k + 1);
        LatticeTester::CopyVect(orbitSeed, lat.orbitSeed, k);
        //      ifm1 = lat.ifm1;
        //      ifr = lat.ifr;
      }
      return *this;
    }


  //===========================================================================

  template<typename Int>
    void MRGComponent<Int>::init (const Int & m0, int k0, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer)
    {
      PolyPE::setM(m0);
      module.init(m0);
      k = k0;
      a.resize(k + 1);
      orbitSeed.resize(k + 1);

      Int m1;
      m1 = getM() - 1;
      ifm1.setNumber (m1);

      if (decom1 == DECOMP_READ)
        ifm1.read (filem1);
      else if (decom1 == DECOMP)
        ifm1.factorize();
      else if (decom1 == DECOMP_WRITE) {
        ifm1.factorize();
        ofstream fout(filem1);

        fout << ifm1.toString();
      }
      ifm1.calcInvFactors();

      Int r;
      r = (power(m0, k) - 1) / (m0 - 1);
      ifr.setNumber(r);

      if (decor == DECOMP_READ)
        ifr.read (filer);
      else if (decor == DECOMP)
        ifr.factorize();
      else if (decor == DECOMP_WRITE) {
        ifr.factorize();
        ofstream fout(filer);
        fout << ifr.toString();
      } else if (decor == DECOMP_PRIME)
        ifr.setStatus (LatticeTester::PRIME);

      ifr.calcInvFactors();
    }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (const Int & m, int k, DecompType decom1,
        const char *filem1, DecompType decor, const char *filer)
    {
      init (m, k, decom1, filem1, decor, filer);
    }


  //===========================================================================

  template<typename Int>
    MRGComponent<Int>::MRGComponent (Modulus<Int> & modu, int k, 
        DecompType decom1, const char *filem1, DecompType decor,
        const char *filer)
    {
      init (modu.m, k, decom1, filem1, decor, filer);

      LatticeTester::PrimeType status = LatticeTester::IntFactor::isPrime (
          modu.m, 100);
      if (status == LatticeTester::PRIME 
          || LatticeTester::PROB_PRIME == status) {
        modu.primeF = true;
      } else {
        cout << " WARNING:  m is NOT prime" << endl;
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
    void MRGComponent<Int>::setA (const MVect & b)
    {
      a = b;
      //  LatticeTester::CopyVect(b, a, k);
    }


  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod (const MVect & a0)
    {
      PolyPE::setM(getM());
      a = a0;
      PolyPE::reverse (a, k, 2);
      PolyPE::setF(a);
      PolyPE pol;
      IntPrimitivity<Int> privfm(ifm1, getM(), 1);
      return pol.isPrimitive(privfm, ifr);
    }


  //===========================================================================

  template<typename Int>
    bool MRGComponent<Int>::maxPeriod23 (const MVect & a0)
    {
      PolyPE::setM(getM());
      a = a0;
      PolyPE::reverse (a, k, 2);
      PolyPE::setF(a);
      PolyPE pol;
      // La condition 1 a déjà été vérifiée dans SeekMain
      return pol.isPrimitive(ifr);
    }


  //===========================================================================

  template<typename Int>
    string MRGComponent<Int>::toString ()
    {
      ostringstream os;
      os << "MRGComponent:";
      Int mm = getM();
      os << "\n   m = " << mm;
      os << "\n   k = " << k;
      os << "\n   a = ";
      string str (os.str ());
      string s2 = LatticeTester::toString(a, k);
      str += s2;
      str += "\n";
      return str;
    }

} // End namespace LatMRG
#endif

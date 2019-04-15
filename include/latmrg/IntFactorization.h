#ifndef INTFACTORIZATION_H
#define INTFACTORIZATION_H
#include <vector>
#include <list>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cassert>

#include "latticetester/Const.h"
#include "latticetester/IntFactor.h"
#include "latticetester/Util.h"

namespace LatMRG {

  /**
   * The class `IntFactorization` implements the decomposition of integers in
   * factors, preferably prime (see class <tt>IntFactor</tt>). It contains
   * functions to factorize an integer in prime factors, to sort and print the
   * list of its factors. `IntFactorization`
   * Integers are factorized by calling the MIRACL software
   * \cite iSCO03a, which uses many different methods in succession to
   * effect the factorization.
   *
   * Given any natural integer \f$n\f$, there is a unique decomposition in
   * prime factors of the form
   * \f[
   *   n = p_1^{\nu_1}p_2^{\nu_2}\cdots p_s^{\nu_s}
   * \f]
   * where \f$p_i\f$ is a prime integer with \f$\nu_i\f$ its multiplicity, and
   * where the factors are sorted in increasing order. In the case of very
   * large integers, it may not be possible to find all the prime factors
   * within a reasonable amount of time. In that case, a similar decomposition
   * to the above may be used with some of the factors composite.
   */
  template <typename Int>
    class IntFactorization {
      public:

        /**
         * Integer \f$x\f$ is the number whose "prime" factors decomposition is
         * kept in this object.
         */
        explicit IntFactorization (const Int & x);

        /**
         * Integer `number` and its prime factors will be read from file
         * `fname` by calling method `read` below. If no argument is given,
         * `number` is initialized to 0.
         */
        explicit IntFactorization (const char *fname = 0);

        /**
         * Destructor.
         */
        ~IntFactorization ();

        /**
         * Copy constructor.
         */
        IntFactorization (const IntFactorization & f);

        /**
         * Assignment operator.
         */
        IntFactorization & operator= (const IntFactorization & f);

        /**
         * Empties the lists of primes factors and set this number to 0.
         */
        void clear ();

        /**
         * Reads the list of (possibly) prime factors of a number from file
         * `f`. The first line contains the number itself. The following lines
         * contain one factor per line: the factor (first field) with its
         * multiplicity (second field) and its status (third field). The status
         * field is written as `P` if the factor is known to be prime, `Q` if
         * the factor is probably prime, `C` if the factor is composite, and
         * `U` if its status is unknown or unimportant. For example, given the
         * number \f$120= 2^3*3*5\f$, then the file must be
         *
         * <center>
         *  <div class="LatSoft-fbox"> <div class="LatSoft-parbox">  120
         * <br>2\  P <br>3  P <br>5  P <br>
         * </div> </div>
         * </center>
         */
        void read (const char *f);

        /**
         * Adds factor \f$p\f$ with multiplicity `mult` and prime status `st`
         * to this object.
         */
        void addFactor (const Int & p, int mult = 1,
            LatticeTester::PrimeType st = LatticeTester::UNKNOWN);

        /**
         * Replaces repeated equal factors in `factorList` by one factor with
         * its multiplicity. Also sorts the factors.
         */
        void unique ();

        /**
         * Tries to find all the prime factors of this number.
         */
        void factorize ();

        /**
         * Checks that the number is equal to the product of its factors.
         * Returns `true` if it is, otherwise `false`.
         */
        bool checkProduct () const;

        /**
         * Given the list of prime factors \f$p\f$ of `number`, fills the list
         * of inverse factors with the values <tt>number</tt>/\f$p\f$.
         */
        void calcInvFactors ();

        /**
         * Returns the value of this number.
         */
        Int getNumber () const { return m_number; }

        /**
         * Returns a non-mutable list of the factors.
         */
        const std::list<LatticeTester::IntFactor<Int>> & getFactorList () const
        { return m_factorList; }

        /**
         * Returns a non-mutable list of the inverse factors.
         */
        const std::vector<Int> & getInvFactorList () const
        { return m_invFactorList; }

        /**
         * Sets the value of this number to \f$x\f$.
         */
        void setNumber (const Int & x) { m_number = x; }

        /**
         * Returns the status of this number.
         */
        LatticeTester::PrimeType getStatus () const { return m_status; }

        /**
         * Sets the status of this number to \f$s\f$.
         */
        void setStatus (LatticeTester::PrimeType s) { m_status = s; }

        /**
         * Returns the list of (possibly) prime factors of this object in the
         * same format as described in method `read` above.
         */
        std::string toString () const;
      private:

        /**
         * The number whose "prime" factor decomposition is kept in this object.
         */
        Int m_number;

        /**
         * The status of this number, i.e. whether it is prime, composite, ...
         */
        LatticeTester::PrimeType m_status;

        /**
         * The list of the "prime" factors in the decomposition of `number`.
         */
        std::list<LatticeTester::IntFactor<Int>> m_factorList;

        /**
         * Given the list of prime factors \f$p\f$ of `number`, `invFactorList`
         * contains all the sorted values <tt>number</tt>/\f$p\f$ (indexing
         * starts at 0). However, one must have called the function
         * `calcInvFactors` beforehand. For example, if `number` = 24, its
         * prime factors decomposition is \f$24 = 2^3\cdot3\f$, and
         * `invfactorList` = \f$[8, 12]\f$.
         */
        std::vector<Int> m_invFactorList;

        /**
         * Nested class used to sort the prime factors of a number in
         * increasing order. Returns `true` if the factor of `f1` is smaller
         * than the factor of `f2`; otherwise returns `false`.
         */
        class CompFactor {
          public:
            bool operator() (const LatticeTester::IntFactor<Int> & f1,
                const LatticeTester::IntFactor<Int> & f2) {
              return f1.getFactor() < f2.getFactor(); }
        };

        /**
         * Sorts the list of factors in increasing order, the smallest factor
         * first.
         */
        void sort () { CompFactor comp; m_factorList.sort (comp); }
    }; // End class IntFactorization

  //IMPLEMENTATIONS

  template<typename Int>
    IntFactorization<Int>::IntFactorization (const char *name):
      m_status(LatticeTester::UNKNOWN)
  {
    if (name != 0)
      read(name);
    else
      m_number = 0;
  }


  //===========================================================================

  template <typename Int>
    IntFactorization<Int>::IntFactorization (const Int & n):
      m_number (n), m_status(LatticeTester::UNKNOWN)
  {}


  //===========================================================================

  template<typename Int>
    IntFactorization<Int>::IntFactorization (const IntFactorization & f)
    : m_number(f.m_number), m_status(f.m_status), m_factorList(f.m_factorList),
    m_invFactorList(f.m_invFactorList)
  {}


  //===========================================================================

  template<typename Int>
    IntFactorization<Int> & IntFactorization<Int>::operator= (
        const IntFactorization & f)
    {
      if (this != &f) {
        m_factorList = f.m_factorList;
        m_invFactorList = f.m_invFactorList;
        m_number = f.m_number;
        m_status = f.m_status;
      }
      return *this;
    }


  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::clear ()
    {
      m_status = LatticeTester::UNKNOWN;
      m_number = 0;
      m_factorList.clear ();
      m_invFactorList.clear ();
    }

  //===========================================================================

  template<typename Int>
    IntFactorization<Int>::~IntFactorization ()
    {
      clear ();
    }


  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::addFactor (const Int & x, int mult,
        LatticeTester::PrimeType st)
    {
      LatticeTester::IntFactor<Int> f (x, mult, st);
      m_factorList.push_back (f);
    }


  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::read (const char *name)
    {
      std::ifstream in (name);

      if (!(in.is_open())) {
        std::string str("IntFactorization::read:   Unable to open input file  ");
        str += name;
        throw std::invalid_argument(str);
      }
      std::string tampon;

      // This should be modified to ignore the entire line because right now we
      // we are limited if the number is too big
      in.ignore (256, '\n'); // drop rest of line

      int vsize = 0;



      while (in >> tampon) {
        Int x;
        int k;
        char c;
        LatticeTester::PrimeType stat;
        NTL::ZZ f;

        f = NTL::to_ZZ (tampon.c_str ());
        NTL::conv (x, f);
        in >> k;
        in >> c;
        switch (c) {

          case 'P':
            stat = LatticeTester::PRIME;
            break;
          case 'Q':
            stat = LatticeTester::PROB_PRIME;
            break;
          case 'C':
            stat = LatticeTester::COMPOSITE;
            break;
          default:
            stat = LatticeTester::UNKNOWN;
        }

        addFactor (x, k, stat);
        ++vsize;
      }

      //unique ();
      assert (checkProduct ());

      m_invFactorList.reserve (vsize);
    }

  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::unique ()
    {
      sort ();
      int j = 1;

      typename std::list<LatticeTester::IntFactor<Int>>::iterator it = m_factorList.begin ();
      typename std::list<LatticeTester::IntFactor<Int>>::iterator it2 = m_factorList.begin ();
      if (it2 != m_factorList.end ())
        ++it2;
      while (it2 != m_factorList.end ()) {
        if (it->getFactor () == it2->getFactor ()) {
          it->setMultiplicity (++j);
          it2 = m_factorList.erase (it2);
        } else {
          j = 1;
          ++it;
          ++it2;
        }
      }
    }


  //===========================================================================

  template<typename Int>
    std::string IntFactorization<Int>::toString () const
    {
      typename std::list<LatticeTester::IntFactor<Int>>::const_iterator it =
        m_factorList.begin ();
      std::ostringstream out;
      out << m_number << std::endl;
      //   out << "its factors:\n";
      while (it != m_factorList.end ()) {
        out << "\t" << (*it).toString () << std::endl;
        ++it;
      }
      out << std::endl;
      /*
         if (!m_invFactorList.empty()) {
         out << "the inverse factors:\n";
         for (unsigned int i = 0; i < m_invFactorList.size(); ++i) {
         out << m_invFactorList[i] << std::endl;
         }
         }
         */
      return out.str ();
    }


  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::factorize ()
    {
#ifdef USE_YAFU
      std::string S("./data/yafu \"factor(");
      std::ostringstream num;
      num << m_number;
      S += num.str() + ")\"";

      // Choose a temporary name for the file
      const char *filename = "temp938573291";

      S += " > ";
      S += filename;
      // factorize and set output to filename
      std::system (S.c_str ());
      // Now read the result file and extract the prime factors from the
      // lines PRIME FACTOR xxx
      std::ifstream in (filename);

      if (!(in.is_open())) {
        std::cerr << "Error:   cannot open file   filename\n";
        exit(8);
      }
      std::string line;
      std::string::size_type pos;
      Int z;
      do {
        getline (in, line, '\n');
        pos = line.find ("factors found");
      } while (pos == std::string::npos);
      while (getline (in, line, '\n')) {
        pos = line.find ("=");
        if (pos != std::string::npos) {
          // Found a prime factor
          S = line.substr (pos + 2);
          NTL::conv(z, S.c_str ());
          addFactor (z, 1, LatticeTester::PRIME);
        }
      }

      unique ();
      remove (filename);
      remove("session.log");
      remove("factor.log");
#else
      std::cout << "IntFactorization: Yafu is not installed in ./data.\n"
        "For more information on how to fix this problem, look at the\n"
        "installation documentation.\n";
      std::cout << "Exiting the program to avoid undefined behavior.";
      exit(1);
#endif
    }

  //===========================================================================

  template<typename Int>
    bool IntFactorization<Int>::checkProduct () const
    {
      Int temp;
      temp = 1;

      typename std::list<LatticeTester::IntFactor<Int>>::const_iterator it =
        m_factorList.begin ();

      while (it != m_factorList.end ()) {
        for (int j = it->getMultiplicity (); j > 0; --j) {
          temp *= it->getFactor ();
        }
        ++it;
      }
      if (temp != m_number) {
        std::cout << "checkProduct:   ERROR -->  " << m_number
          << " != " << temp << std::endl;
      }
      return (temp == m_number);
    }

  //===========================================================================

  template<typename Int>
    void IntFactorization<Int>::calcInvFactors ()
    {
      //   sort ();
      unsigned int j = 0;

      for (auto it = m_factorList.rbegin(); it != m_factorList.rend(); it++) {
        if (it->getFactor() == Int(1)) continue;
        ++j;
        if (m_invFactorList.capacity () < j) {
          m_invFactorList.clear ();
          m_invFactorList.reserve (j + 10);
        }
        m_invFactorList.push_back (m_number / it->getFactor ());
      }
    }

  extern template class IntFactorization<std::int64_t>;
  extern template class IntFactorization<NTL::ZZ>;

}
#endif

#ifndef LATMRG_INTFACTORIZATION_H
#define LATMRG_INTFACTORIZATION_H
#define USE_YAFU

#include <vector>
#include <list>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cassert>

// #include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latmrg/IntFactor.h"

namespace LatMRG {

/**
 * Represents the factorization of an arbitrary positive integer, usually into prime factors,
 * but not always.  The factors are `IntFactor` objects.
 * There are also functions to sort and print the list of factors.
 * The factorization is performed by the `MIRACL` software \cite iSCO03a,
 * which uses different methods in succession to perform the factorization.
 *
 * Recall that for any natural integer \f$n\f$, there is a unique decomposition in
 * prime factors of the form
 * \f[
 *   n = p_1^{\nu_1} p_2^{\nu_2} \cdots p_s^{\nu_s}
 * \f]
 * where \f$p_i\f$ is a prime factor with multiplicity \f$\nu_i\f$ and
 * the factors are sorted in increasing order.  For very
 * large integers, it may not be possible to find all the prime factors
 * within a reasonable amount of time. In that case, a similar decomposition
 * to the above may be used with some of the factors composite.
 */
template<typename Int>
class IntFactorization {
public:

   /**
    * Constructs a factorization for the integer \f$x\f$.
    */
   explicit IntFactorization(const Int &x);

   /**
    * Constructs a factorization object by reading the factorization from file `fname`
    * using the `read` function defined below. If no file name is given, the integer is
    * initialized to 0.
    */
   explicit IntFactorization(const char *fname = 0);

   /**
    * Destructor.
    */
   ~IntFactorization();

   /**
    * Copy constructor.
    */
   IntFactorization(const IntFactorization &f);

   /**
    * Assignment operator.
    */
   IntFactorization& operator=(const IntFactorization &f);

   /**
    * Empties the list of primes factors and set this number to 0.
    */
   void clear();

   /**
    * Reads a factorization by reading the list of (possibly prime) factors
    * of an integer from file `f`, in the following format.
    * The first line contains the integer itself. The following lines
    * contain one factor per line: the factor (first field) with its
    * multiplicity (second field), and its status (third field). The status
    * field is written as `P` if the factor is known to be prime, `Q` if
    * the factor is probably prime, `C` if the factor is composite, and
    * `U` if its status is unknown (or unimportant). For example, for the
    * number \f$x = 120 = 2^3*3*5\f$, the file should contain
    *
    * <center>
    *  <div class="LatSoft-fbox"> <div class="LatSoft-parbox">
    *  120 <br>2 1 P <br>3 2 P <br>5 1 P <br>
    * </div> </div>
    * </center>
    */
   void read(const char *f);

   /**
    * Adds the factor `p` with multiplicity `mult` and prime status `status`
    * to this object.
    */
   void addFactor(const Int &p, int64_t mult = 1, PrimeType status = UNKNOWN);

   /**
    * Replaces repeated equal factors in the factor list of this object by one factor with
    * its multiplicity, and sorts the factors in increasing order.
    */
   void makeUnique();

   /**
    * Tries to find all the prime factors of this integer.
    * The current implementation uses Yafu and works only if Yafu is installed.
    */
   void factorize();

   /**
    * Checks that the integer in this object is equal to the product of its factors.
    * Returns `true` if it is, otherwise `false`.
    */
   bool checkProduct() const;

   /**
    * Given the list of prime factors \f$p\f$ of the integer \f$x\f$ in this object,
    * this function computes the internal list of the inverse factors \f$x/p\f$,
    * sorted in increasing order.
    * For example, if \f$x = 24\f$, the prime factorization is \f$24 = 2^3 \cdot 3\f$
    * and the list on inverse factors is \f$[8, 12]\f$.
    */
   void calcInvFactors();

   /**
    * Returns the main integer in this object.
    */
   Int getNumber() const {
      return m_number;
   }

   /**
    * Returns the list of factors.
    */
   const std::list<IntFactor<Int>>& getFactorList() const {
      return m_factorList;
   }

   /**
    * Returns the list of the inverse factors, assuming that they have been computed.
    */
   const std::vector<Int>& getInvFactorList() const {
      return m_invFactorList;
   }

   /**
    * Sets the value of this integer to \f$x\f$.
    */
   void setNumber(const Int &x) {
      m_number = x;
   }

   /**
    * Returns the status of this number.
    */
   PrimeType getStatus() const {
      return m_status;
   }

   /**
    * Sets the status of this number to `status`.
    */
   void setStatus(PrimeType status) {
      m_status = status;
   }

   /**
    * Returns the list of (possibly prime) factors of this object as a string in the
    * same format as described in method `read` above.
    */
   std::string toString() const;

private:

   /**
    * The number whose factor decomposition is kept in this object.
    */
   Int m_number;

   /**
    * The status of this number, i.e. whether it is prime, composite, ...
    */
   PrimeType m_status;

   /**
    * The list of factors in the decomposition of number.
    */
   std::list<IntFactor<Int>> m_factorList;

   /**
    * Given the list of factors \f$p\f$ of this integer \f$x\f$, `m_invFactorList`
    * contains all the sorted values \f$x/p\f$ (indexing starts at 0),
    * if the function `calcInvFactors` has been called beforehand.
    */
   std::vector<Int> m_invFactorList;

   /**
    * Nested class used to sort the prime factors of an integer in
    * increasing order. Returns `true` if the factor `f1` is smaller
    * than the factor `f2`; otherwise returns `false`.
    */
   class CompareFactors {
   public:
      bool operator()(const IntFactor<Int> &f1, const IntFactor<Int> &f2) {
         return f1.getFactor() < f2.getFactor();
      }
   };

   /**
    * Sorts the list of factors in increasing order.
    */
   void sort() {
      CompareFactors comp;
      m_factorList.sort(comp);
   }

};
// End class IntFactorization

// ==============================================================================

// IMPLEMENTATION

template<typename Int>
IntFactorization<Int>::IntFactorization(const char *name) :
      m_status(UNKNOWN) {
   if (name != 0) read(name);
   else m_number = 0;
}

//===========================================================================

template<typename Int>
IntFactorization<Int>::IntFactorization(const Int &n) :
      m_number(n), m_status(UNKNOWN) {
}

//===========================================================================

template<typename Int>
IntFactorization<Int>::IntFactorization(const IntFactorization &f) :
      m_number(f.m_number), m_status(f.m_status), m_factorList(f.m_factorList), m_invFactorList(
            f.m_invFactorList) {
}

//===========================================================================

template<typename Int>
IntFactorization<Int>& IntFactorization<Int>::operator=(const IntFactorization &f) {
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
void IntFactorization<Int>::clear() {
   m_status = UNKNOWN;
   m_number = 0;
   m_factorList.clear();
   m_invFactorList.clear();
}

//===========================================================================

template<typename Int>
IntFactorization<Int>::~IntFactorization() {
   clear();
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::addFactor(const Int &x, int64_t mult, PrimeType st) {
   IntFactor<Int> f(x, mult, st);
   m_factorList.push_back(f);
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::read(const char *name) {
   std::ifstream in(name);
   if (!(in.is_open())) {
      std::string str("IntFactorization::read:   Unable to open input file  ");
      str += name;
      throw std::invalid_argument(str);
   }
   std::string tampon;
   // This should be modified to ignore the entire line because right now
   // we are limited if the number is too big.                             *********
   in.ignore(256, '\n'); // drop rest of line
   int64_t vsize = 0;
   while (in >> tampon) {
      Int x;
      int64_t k;
      char c;
      PrimeType status;
      NTL::ZZ f;
      f = NTL::to_ZZ(tampon.c_str());
      NTL::conv(x, f);
      in >> k;
      in >> c;
      switch (c) {
      case 'P':
         status = PRIME;
         break;
      case 'Q':
         status = PROB_PRIME;
         break;
      case 'C':
         status = COMPOSITE;
         break;
      default:
         status = UNKNOWN;
      }
      addFactor(x, k, status);
      ++vsize;
   }

   // makeUnique ();
   assert(checkProduct());
   m_invFactorList.reserve(vsize);
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::makeUnique() {
   sort();
   int64_t j = 1;
   typename std::list<IntFactor<Int>>::iterator it = m_factorList.begin();
   typename std::list<IntFactor<Int>>::iterator it2 = m_factorList.begin();
   if (it2 != m_factorList.end()) ++it2;
   while (it2 != m_factorList.end()) {
      if (it->getFactor() == it2->getFactor()) {
         it->setMultiplicity(++j);
         it2 = m_factorList.erase(it2);
      } else {
         j = 1;
         ++it;
         ++it2;
      }
   }
}

//===========================================================================

template<typename Int>
std::string IntFactorization<Int>::toString() const {
   typename std::list<IntFactor<Int>>::const_iterator it = m_factorList.begin();
   std::ostringstream out;
   out << m_number << ":\n";
   //   out << "its factors:\n";
   while (it != m_factorList.end()) {
      out << (*it).toString() << std::endl;
      ++it;
   }
   out << std::endl;
   /*
    if (!m_invFactorList.empty()) {
    out << "the inverse factors:\n";
    for (unsigned int64_t i = 0; i < m_invFactorList.size(); ++i) {
    out << m_invFactorList[i] << std::endl;
    }
    }
    */
   return out.str();
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::factorize() {
#ifdef USE_YAFU
      //std::string S("./data/yafu \"factor(");
      std::string S("./data/yafu -s ");
      std::ostringstream num;
      num << m_number;
      //S += num.str() + ")\"";
      
      S += num.str();
      
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
      //std::string::size_type pos;
      Int z;
      
      
      while (getline (in, line)) {
        S = line;
        NTL::conv(z, S.c_str ());
        addFactor (z, 1, PRIME);
      }
      
      makeUnique();
      remove (filename);
      
      /*
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
          addFactor (z, 1, PRIME);
        }
      }
      
      unique ();
      remove (filename);
      remove("session.log");
      remove("factor.log");
      remove("siqs.dat");
      */
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
bool IntFactorization<Int>::checkProduct() const {
   Int temp;
   temp = 1;
   typename std::list<IntFactor<Int>>::const_iterator it = m_factorList.begin();
   while (it != m_factorList.end()) {
      for (int64_t j = it->getMultiplicity(); j > 0; --j) {
         temp *= it->getFactor();
      }
      ++it;
   }
   if (temp != m_number) {
      std::cout << "checkProduct: ERROR -->  " << m_number << " != " << temp << std::endl;
   }
   return (temp == m_number);
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::calcInvFactors() {
   //   sort ();
   // unsigned
   uint64_t j = 0;
   for (auto it = m_factorList.rbegin(); it != m_factorList.rend(); it++) {
      if (it->getFactor() == Int(1)) continue;
      ++j;
      if (m_invFactorList.capacity() < j) {
         m_invFactorList.clear();
         m_invFactorList.reserve(j + 10);
      }
      m_invFactorList.push_back(m_number / it->getFactor());
   }
}

template class IntFactorization<std::int64_t> ;
template class IntFactorization<NTL::ZZ> ;

}
#endif

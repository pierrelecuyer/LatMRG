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
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"

namespace LatMRG {

/**
 * Represents the factorization of an arbitrary positive integer, usually into prime factors,
 * but not always.  The factors are `IntFactor` objects.
 * There are also functions to sort and print the list of factors.
 * The factorization is performed by the `yafu` executable program, which must be installed
 * in a directory that belongs to the `PATH` for the function to work.
 * This program is in `latmrg/data/yafu` in the git distribution.
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
 *
 * The functions that perform a decomposition and return a `PrimeType` will return
 * `PRIME` if all factors are known to be prime,
 * `PROB_PRIME` if they are all prime or probably prime,
 * `COMPOSITE` if at least one factor is known to be composite and could not be decomposed,
 * and `UNKNOWN` in other cases.
 * This returned type is called the *status of the factorization*, and is also returned by
 * the function `getFactStatus`. A `PROB_PRIME` status means that we have a correct factorization
 * into prime numbers with high probability, not that the number \f$n\f$ itself is prime with
 * high probability.  For the latter, there is the function `getNumberStatus` which returns
 * `PRIME` if \f$n\f$ is known to be prime, `PROB_PRIME` if it is probably prime, etc.
 */
template<typename Int>
class IntFactorization {

public:

   /**
    * Constructs a factorization for the natural integer \f$n\f$.
    */
   explicit IntFactorization(const Int &n);

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
    * of an integer from file `f` whose name is given as a character string, in the following format.
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
   void makeUniqueAndSort();

   /**
    * Tries to find all the prime factors of this integer and stores all the factors
    * in increasing order in a private list.  Then calls `makeUniqueAndSort`.
    * If the number is prime, there is a single factor in the list.
    * This function works only if Yafu is installed.
    * Returns the status of the factorization, which is
    * `PRIME` if all factors are known to be prime,
    * `PROB_PRIME` if they are all prime or probably prime,
    * `COMPOSITE` if at least one factor is known to be composite and could not be decomposed,
    * and `UNKNOWN` in other cases.
    */
   PrimeType factorize();

   /**
    * Given the list of prime factors \f$p\f$ of the integer \f$n\f$ in this object,
    * this function computes the internal list of the inverse factors \f$n/p\f$,
    * sorted in increasing order.
    * For example, if \f$n = 24\f$, the prime factorization is \f$24 = 2^3 \cdot 3\f$
    * and the list on inverse factors is \f$[24/3, 24/2] = [8, 12]\f$.
    */
   void calcInvFactors();

   /**
    * Makes a decomposition of the current integer according to the selected value of `decomp`.
    * The type `DecompType` is defined in `EnumTypes`.  The choices are
    * DECOMP, DECOMP_WRITE, DECOMP_READ, DECOMP_PRIME, NO_DECOMP'.
    * Unless this value is `NO_DECOMP`, the prime factors as well as the inverse factors
    * are stored in the current object.
    * The string `filename` is the name of the file where the factors are read or written
    * when `decomp` is <tt>DECOMP_WRITE</tt> or <tt>DECOMP_READ</tt>.
    * This file must be accessible by the program.
    * Return is like in `factorize`.
    */
   PrimeType decompToFactorsInv(DecompType decomp, const char *file = NULL);

   /**
    * Same as calling `factorize` and `calcInvFactors` in a single function call.
    * Return is like in `factorize`.
    */
   PrimeType factorizePlus();

   /**
    * Checks that the integer in this object is equal to the product of its factors.
    */
   bool checkProduct() const;

   /**
    * Returns the main integer \f$n\f$ in this object.
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
    * Returns the vector of the inverse factors, assuming that they have been computed.
    */
   const std::vector<Int>& getInvFactorList() const {
      return m_invFactorList;
   }

   /**
    * Sets to \f$n\f$ the value of the integer to be factored, and clears the list of factors.
    */
   void setNumber(const Int &n) {
      clear();
      m_number = n;
   }

   /**
    * Returns the status of this number.
    */
   PrimeType getNumberStatus() const {
      return m_numberStatus;
   }

   /**
    * Sets the status of this number to `status`.
    */
   void setNumberStatus(PrimeType status) {
      m_numberStatus = status;
   }

   /**
    * Returns the status of this factorization.
    */
   PrimeType getFactStatus() const {
      return m_factStatus;
   }

   /**
    * Sets the status of this factorization to `status`.
    */
   void setFactStatus(PrimeType status) {
      m_factStatus = status;
   }

   /**
    * Returns the list of (possibly prime) factors of this object as a string in the
    * same format as described in method `read` above.
    */
   std::string toString() const;

private:

   /**
    * The integer whose factor decomposition is kept in this object.
    */
   Int m_number;

   /**
    * The status of the natural integer \f$n\f$.
    */
   PrimeType m_numberStatus;

   /**
    * The status of the current decomposition of this number, i.e. the `PrimeType`
    * that was returned by the function that decomposed it for the last time.
    */
   PrimeType m_factStatus;

   /**
    * The list of factors in the decomposition of number.
    */
   std::list<IntFactor<Int>> m_factorList;

   /**
    * Given the list of factors \f$p\f$ of this integer \f$x\f$, the `m_invFactorList`
    * vector contains all the sorted values \f$x/p\f$ (indexing starts at 0),
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
   void sortFactors() {
      CompareFactors comp;
      m_factorList.sort(comp);
   }

};
// End class IntFactorization

// ==============================================================================
// IMPLEMENTATION

template<typename Int>
IntFactorization<Int>::IntFactorization(const char *name) :
      m_factStatus(UNKNOWN) {
   if (name != 0) read(name);
   else m_number = 0;
}

//===========================================================================

template<typename Int>
IntFactorization<Int>::IntFactorization(const Int &n) :
      m_number(n), m_factStatus(UNKNOWN) {
}

//===========================================================================

template<typename Int>
IntFactorization<Int>::IntFactorization(const IntFactorization &f) :
      m_number(f.m_number), m_numberStatus(f.m_numberStatus), m_factStatus(f.m_factStatus), m_factorList(
            f.m_factorList), m_invFactorList(f.m_invFactorList) {
}

//===========================================================================

template<typename Int>
IntFactorization<Int>& IntFactorization<Int>::operator=(const IntFactorization &f) {
   if (this != &f) {
      m_factorList = f.m_factorList;
      m_invFactorList = f.m_invFactorList;
      m_number = f.m_number;
      m_numberStatus = f.m_numberStatus;
      m_factStatus = f.m_factStatus;
   }
   return *this;
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::clear() {
   m_numberStatus = UNKNOWN;
   m_factStatus = UNKNOWN;
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
void IntFactorization<Int>::read(const char *name) {
   std::cout << "Reading file of factors " << name << "\n";
   std::ifstream in(name);
   if (!(in.is_open())) {
      std::string str("IntFactorization::read: Unable to open input file  ");
      str += name;
      throw std::invalid_argument(str);
   }
   std::cout << "Declaring string `tampon` \n";
   std::string tampon;
   // This should be modified to ignore the entire line because right now
   // we are limited if the number is too big.   ?????
   // in.ignore(256, '\n'); // drop rest of line
   int64_t vsize = 0;
   in >> m_number;
   std::cout << "Number: " << m_number << "\n";
   m_factStatus = PRIME;
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
      m_factStatus = max(status, m_factStatus);
      m_numberStatus = UNKNOWN;
      ++vsize;
   }
   // makeUniqueAndSort ();
   assert(checkProduct());
   m_invFactorList.reserve(vsize);
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::addFactor(const Int &x, int64_t mult, PrimeType st) {
   IntFactor<Int> f(x, mult, st);
   m_factorList.push_back(f);
   m_factStatus = max(st, m_factStatus);
   if (m_factorList.size() > 1) m_numberStatus = COMPOSITE;
}

//===========================================================================

template<typename Int>
PrimeType IntFactorization<Int>::decompToFactorsInv(DecompType decomp, const char *filename) {
   if (decomp == NO_DECOMP) return m_factStatus = UNKNOWN;
   if (decomp == DECOMP_READ) {
      read(filename);
      m_factStatus = PRIME;
   } else m_factStatus = factorize();
   if ((decomp == DECOMP_WRITE) && (m_factStatus <= 1)) {
      std::ofstream fout(filename);
      fout << toString();
   }
   if (m_factStatus <= 1) calcInvFactors();
   return m_factStatus;
}

//===========================================================================

template<typename Int>
PrimeType IntFactorization<Int>::factorize() {
#if defined(USE_YAFU)
   // std::cout << "  Start factorize with Yafu " << "\n";
   std::ostringstream num;
   num << m_number;
 
   std::string command = "yafu \"factor(" + num.str() + ")\"";
   FILE* pipe = popen(command.c_str(), "r");
   if (!pipe) throw std::runtime_error("Failed to start YAFU");

   std::string line;
   std::string::size_type pos;
   Int z;
   m_factStatus = PROB_PRIME;
   bool foundMarker = false;
   char buffer[4096];
   
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
      line = buffer;
      // Delete "\n" if this is at the end of the string
      if (!line.empty() && line.back() == '\n') line.pop_back(); 

      // Check if we are in the relevant part of the yafu output
      if (!foundMarker) {
            if (line.find("***factors found***") != std::string::npos) {
                foundMarker = true; 
            } else {
                continue; 
            }
      }
      // Only read out output lines which start with P (for prime) or C (for composite number)
      if (!line.empty() && (line[0] == 'P' || line[0] == 'C')) {         
          pos = line.find("= ");
          if (pos == std::string::npos) {
             throw std::runtime_error("Unexpected YAFU output format");
          }
          NTL::conv(z, line.substr(pos + 2).c_str());
          if (z != 0 && line[0] == 'P') addFactor(z, 1, PROB_PRIME);
          if (z != 0 && line[0] == 'C') {
             addFactor(z, 1, COMPOSITE);
             m_factStatus = COMPOSITE;
          }
      }
      // The last relevant line of the yafu output starts with 'ans = '
      if (line.rfind("ans = ", 0) == 0) break; 
   }
   // Check if the last line is 'ans = 1' which indicates that yafu found a correct decomposition.
   pos = line.find("= ");
   if (pos == std::string::npos || line.substr(pos + 2) != "1") {
      std::cout << "yafu output indicates incomplete prime decomposition \n \n";
   }
   // std::cout << "  the factors are now in m_factorList " << "\n";
   if (m_factorList.size() == 1 && m_factStatus != COMPOSITE) m_numberStatus = PROB_PRIME;
   else m_numberStatus = COMPOSITE;
   makeUniqueAndSort();
   if (pclose(pipe) != 0) {
     std::cerr << "Warning: YAFU returned an error code\n \n";
   }
   // std::cout << "factorize: m_factStatus = " << m_factStatus << "\n";
   return m_factStatus;
#else
   std::cout << "IntFactorization: Yafu is not installed or not accessible.\n";
   std::cout << "Exiting the program to avoid undefined behavior.";
   exit(1);
#endif
}

//===========================================================================

template<typename Int>
PrimeType IntFactorization<Int>::factorizePlus() {
   m_factStatus = this->factorize();
   // std::cout << "factorizePlus: m_factStatus = " << m_factStatus << "\n";
   if (m_factStatus <= 1) calcInvFactors();
   // std::cout << "factorizePlus: after calcInvFactors, size = " << m_invFactorList.size() << "\n";
   return m_factStatus;
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::calcInvFactors() {
   //   sortFactors ();
   if (m_invFactorList.capacity() < m_factorList.size()) {
      m_invFactorList.clear();
      m_invFactorList.reserve(m_factorList.size());
   }
   uint64_t j = 0;
   for (auto it = m_factorList.rbegin(); it != m_factorList.rend(); it++) {
      if (it->getFactor() == Int(1)) continue;
      ++j;
      m_invFactorList.push_back(m_number / it->getFactor());
   }
   // std::cout << "Inside calcInvFactors, size = " << m_invFactorList.size() << "\n";
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
      std::cout << "checkProduct: ERROR --> m_number = " << m_number << " != " << temp << std::endl;
   }
   return (temp == m_number);
}

//===========================================================================

template<typename Int>
void IntFactorization<Int>::makeUniqueAndSort() {
   sortFactors();
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
   out << m_number << "\n";
   while (it != m_factorList.end()) {
      out << (*it).toString() << std::endl;
      ++it;
   }
   // out << std::endl;
   if (!m_invFactorList.empty()) {
      out << "the inverse factors:\n";
      for (uint64_t i = 0; i < m_invFactorList.size();  ++i)
            out << m_invFactorList[i] << std::endl;
   }
   return out.str();
}

// template class IntFactorization<std::int64_t> ;
// template class IntFactorization<NTL::ZZ> ;

   }
#endif

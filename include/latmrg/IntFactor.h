// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATMRG_INTFACTOR_H
#define LATMRG_INTFACTOR_H

#include <string>
#include <iomanip>
#include <cstdint>
#include <sstream>
#include <cassert>
#include <NTL/ZZ.h>

// #include "latticetester/EnumTypes.h"
#include "latmrg/EnumTypes.h"

namespace LatMRG {

static constexpr uint64_t NB_PRIMES = 6543;
const std::array<uint64_t, NB_PRIMES> PRIMES_ARRAY = {
#include "../data/primes.dat"
};

/**
 * An object of this class represents a factor in the decomposition of a positive integer.
 * It is usually a prime factor, but not always.
 * The file also contains basic static functions to test whether an integer is prime,
 * probably prime, or composite. The primality is tested by first checking a list of
 * small prime factors, and then invoking the `NTL::ProbPrime` function.
 */
template<typename Int> class IntFactor {

public:

   /**
    * Constructs a factor \f$x\f$ of multiplicity `mult` with given `PrimeType status`.
    */
   IntFactor(const Int &x, int64_t mult = 1, PrimeType status = UNKNOWN) :
         m_factor(x), m_multiplicity(mult), m_status(status) {
   }

   /**
    * Returns the numerical value of this factor.
    */
   Int getFactor() const {
      return m_factor;
   }

   /**
    * Sets the value of this factor to \f$x\f$.
    */
   void setFactor(const Int &x) {
      m_factor = x;
   }

   /**
    * Returns the multiplicity of this factor.
    */
   int64_t getMultiplicity() const {
      return m_multiplicity;
   }

   /**
    * Sets the multiplicity of this factor to \f$m\f$.
    */
   void setMultiplicity(int64_t m) {
      m_multiplicity = m;
   }

   /**
    * Returns the `PrimeType` of this factor.
    */
   PrimeType getStatus() const {
      return m_status;
   }

   /**
    * Sets the `PrimeType` of this factor to \f$s\f$.
    */
   void setStatus(PrimeType s) {
      m_status = s;
   }

   /**
    * Same as `typePrime (y, numtrials)` with the current object in place of `y`.
    */
   PrimeType typePrime(std::int64_t numtrials = 100);

   /**
    * Same as `isPrime (y, numtrials)` with the current object in place of `y`.
    */
   bool isPrime(std::int64_t numtrials = 100);

   /**
    * Tests whether \f$y\f$ is prime or not, via `NTL::ProbPrime`.
    * The latter first tests whether \f$y\f$ is divisible by small primes.
    * If no factor is found, it applies the Miller-Rabin probability
    * test with `numtrials` trials.  Returns a `PrimeType`, NOT a boolean.
    */
   static PrimeType typePrime(const Int &y, std::int64_t numtrials = 100);

   /**
    * Tests whether \f$y\f$ is prime or not, via `NTL::ProbPrime`.
    * If no small factor is found, the Miller-Rabin probability
    * test from NTL is applied with `numtrials` trials.
    * Returns `true` iff \f$y\f$ is prime or probably prime.
    */
   static bool isPrime(const Int &y, std::int64_t numtrials = 100);

   /**
    * Same as `isPrime`.
    */
   static bool isProbPrime(const Int &y, std::int64_t numtrials = 100);

   /**
    * Similar to `isPrime`, but tests whether \f$y\f$ is a safe prime or not.
    * This means that  \f$y\f$ is prime and  \f$(y-1)/2\f$ is also prime.
    * In this case, we also say that \f$(y-1)/2\f$ is a Sophie Germain prime number.
    */
   static bool isSafePrime(const Int &y, std::int64_t numtrials= 100);

   /**
    * Transforms the status `status` to a string and returns it.
    */
   static std::string toString(PrimeType status);

   /**
    * Returns a string that represents this object.
    */
   std::string toString() const;



   //===========================================================================

private:

   /**
    * This factor as an `Int`.
    */
   Int m_factor;

   /**
    * The multiplicity of this factor.
    */
   int64_t m_multiplicity;

   /**
    * The status of this factor.
    */
   PrimeType m_status;

};
// class IntFactor

//===========================================================================

template<typename Int>
std::string IntFactor<Int>::toString() const {
   char c;
   switch (m_status) {
   case PRIME:
      c = 'P';
      break;
   case PROB_PRIME:
      c = 'Q';
      break;
   case COMPOSITE:
      c = 'C';
      break;
   default:
      c = 'U';
      break;
   }
   std::ostringstream sortie;
   sortie << m_factor << std::setw(10) << m_multiplicity << std::setw(10) << c;
   return sortie.str();
}

//===========================================================================

template<typename Int>
inline std::string IntFactor<Int>::toString(PrimeType status) {
   return toStringPrimeType(status);
}

//===========================================================================

template<typename Int>
PrimeType IntFactor<Int>::typePrime(std::int64_t numtrials) {
   if (isPrime(m_factor, numtrials)) {
      return m_status = PROB_PRIME;
   }
   else {
      return m_status = COMPOSITE;
   }
}

//===========================================================================

template<typename Int>
bool IntFactor<Int>::isPrime(std::int64_t numtrials) {
   if (isPrime(m_factor, numtrials)) {
      m_status = PROB_PRIME;
      return true;
   }
   else {
      m_status = COMPOSITE;
      return false;
   }
}

//===========================================================================

template<typename Int>
PrimeType IntFactor<Int>::typePrime(const Int &y, std::int64_t numtrials) {
   if (isPrime(y, numtrials))
      return PROB_PRIME;
   else
      return COMPOSITE;
}

//===========================================================================

template<typename Int>
inline bool IntFactor<Int>::isPrime(const Int &y, std::int64_t numtrials) {
   return NTL::ProbPrime(y, numtrials);    // Returns 1 if prob. prime.
}

//===========================================================================

template<typename Int>
inline bool IntFactor<Int>::isProbPrime(const Int &y, std::int64_t numtrials) {
   return NTL::ProbPrime(y, numtrials);    // Returns 1 if prob. prime.
}

//===========================================================================

template<typename Int>
inline bool IntFactor<Int>::isSafePrime(const Int &y, std::int64_t numtrials) {
   if (!isPrime(y, numtrials)) return false;
   return isPrime((y-Int(1))/Int(2), numtrials);
}

// template class IntFactor<NTL::ZZ> ;
// template class IntFactor<std::int64_t> ;

} // namespace LatMRG

#endif

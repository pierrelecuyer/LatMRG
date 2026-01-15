// This file is part of LatMRG.
//
// Copyright (C) 2012-2024  The LatMRG authors, under the supervision
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

#ifndef LATMRG_PrimesFinder_H
#define LATMRG_PrimesFinder_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <ctime>

#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
#include "latmrg/IntFactor.h"
#include "latmrg/IntFactorization.h"

/**
 * \file latmrg/PrimesFinder.h
 *
 * This file provides static functions to search for integers \f$m\f$ that are prime and
 * for which the integer \f$r = (m^k-1)/(m-1)\f$ is also prime for a given
 * \f$k\f$, and possibly for which \f$(m-1)/2\f$ is also prime.
 * For \f$k=1\f$, we have \f$r=1\f$, so the first condition is removed.
 * These functions use solely Miller-Rabin probabilistic primality tests with `numtrials` trials.
 * If `factomm1` is `true`, the factorization of \f$m-1\f$ is also computed and returned.
 * Each function writes its results on the given stream `fout`.
 *
 * We could possibly use deterministic primality tests when \f$m\f$ and \f$k\f$ are small, because it is possible
 * to do everything this class does with efficient deterministic tests for integers \f$ < 2^{64}\f$;
 * see https://math.stackexchange.com/questions/2481148/primality-testing-for-64-bit-numbers.
 * However \f$r\f$ is typically much larger than \f$2^{64}\f$ in actual cases.
 */

using namespace LatticeTester;

namespace LatMRG {

/**
 * Finds and prints the \f$s\f$ largest prime integers \f$m < 2^e\f$.
 */
template<typename Int>
static void findPrime(int64_t e, int64_t s, bool factomm1, std::ostream &fout,
      const int64_t numtrials = 200);
/**
 * Finds the \f$s\f$ largest prime integers \f$m<2^e\f$ for which \f$r = (m^k-1)/(m-1)\f$ is also prime.
 * If `safe` is `true`, \f$(m-1)/2\f$ is also required to be prime.
 */
template<typename Int>
static void findPrime(int64_t k, int64_t e, int64_t s, bool safe, bool factomm1, std::ostream &fout,
      const int64_t numtrials = 200);

/**
 * Finds all integers \f$m\f$, in \f$2^e + c_1 \le m \le 2^e + c_2\f$,
 * such that \f$m\f$ and \f$r = (m^k-1)/(m-1)\f$ are prime. If `safe`
 * is `true`, \f$(m-1)/2\f$ is also required to be prime.
 */
template<typename Int>
static void findPrime(int64_t k, int64_t e, int64_t c1, int64_t c2, bool safe, bool factomm1,
      std::ostream &fout, const int64_t numtrials = 200);

/**
 * This is the general purpose function, called by the other ones.
 * It searches for the `s` largest prime integers between `S1` and `S2`.
 */
template<typename Int>
static void findPrime(int64_t k, int64_t e, int64_t s, const Int &S1, const Int &S2, bool safe,
      bool factomm1, std::ostream &fout, const int64_t numtrials = 200);

/**
 * Writes the search parameters to the stream `fout`.
 */
static void writeHeader(int64_t k, int64_t e, int64_t c1, int64_t c2, bool safe, bool factomm1,
      std::ostream &fout, const int64_t numtrials = 200);

//============================================================================
// Implementation

// This is used in the next function.
template<typename Int>
static void nextM(Int &m) {
   m -= 2;
   if (0 == m % 5) m -= 2;
}
;

// This implements the general case.
template<typename Int>
void findPrime(int64_t k, int64_t e, int64_t s, const Int &S1, const Int &S2, bool safe,
      bool factomm1, std::ostream &fout, const int64_t numtrials) {
   Chrono timer;
   timer.init();
   Int m;
   if (NTL::IsOdd(S2)) m = S2;
   else m = S2 - Int(1);
   int64_t i = 0;
   bool rprime;
   while (i < s && m >= S1) {
      if (IntFactor<Int>::isPrime(m, numtrials)) {
         if (factomm1) fout << "----------------\n";
         Int m1 = m - Int(1);
         if (safe) {
            if (1 == m % 4) {
               nextM(m);
               continue;
            }
            Int m1s2 = m1 / Int(2);
            if (!IntFactor<Int>::isPrime(m1s2, numtrials)) {
               nextM(m);
               continue;
            }
         }
         Int r;
         NTL::set(r);
         if (k == 1) rprime = true;
         if (k > 1) {
            r = NTL::power(m, k);
            --r;
            r = r / m1;
            rprime = IntFactor<Int>::isPrime(r, numtrials);
         }
         if (rprime) {
            i++;
            fout << "m = " << m;
            Int Sdiff = m - (Int(1) << e);
            fout << "  = 2^" << e;
            if (Sdiff >= 0) {
               fout << " + ";
               fout << Sdiff << "\n";
            } else {
               fout << " - ";
               fout << (-Sdiff) << "\n";
            }
            if (factomm1) {
               IntFactorization<Int> ifac;
               ifac.clear();
               ifac.setNumber(m1);
               ifac.factorize();
               fout << "Factors of m - 1 = ";
               fout << ifac.toString();
            }
         }
      }
      nextM(m);
   }
   fout << "\nCPU time: " << timer.toString() << "\n\n";
}

//===========================================================================

template<typename Int>
void findPrime(int64_t e, int64_t s, bool factomm1, std::ostream &fout, const int64_t numtrials) {
   Int Sm1, Sm2;
   Sm2 = (Int(1) << e) - 1;
   Sm1 = 2;
   writeHeader(1, e, INT_MAX, INT_MAX, false, factomm1, fout);
   findPrime(1, e, s, Sm1, Sm2, false, factomm1, fout, numtrials);
}

//===========================================================================

template<typename Int>
void findPrime(int64_t k, int64_t e, int64_t s, bool safe, bool factomm1, std::ostream &fout,
      const int64_t numtrials) {
   Int Sm1, Sm2;
   Sm2 = (Int(1) << e) - 1;
   Sm1 = 2;
   writeHeader(k, e, INT_MAX, INT_MAX, safe, factomm1, fout);
   findPrime(k, e, s, Sm1, Sm2, safe, factomm1, fout, numtrials);
}

//===========================================================================

template<typename Int>
void findPrime(int64_t k, int64_t e, int64_t c1, int64_t c2, bool safe, bool factomm1,
      std::ostream &fout, const int64_t numtrials) {
   Int Sm1, Sm2;
   Sm1 = (Int(1) << e) + c1;
   Sm2 = (Int(1) << e) + c2;
   writeHeader(k, e, c1, c2, safe, factomm1, fout);
   findPrime(k, e, INT_MAX, Sm1, Sm2, safe, factomm1, fout, numtrials);
}

//===========================================================================

void writeHeader(int64_t k, int64_t e, int64_t c1, int64_t c2, bool safe, bool factomm1,
      std::ostream &fout, const int64_t numtrials) {
   fout << "//=============================================================\n";
   fout << "Largest values of m such that m";
   if (safe) fout << ", (m-1)/2,";
   if (k > 1) fout << " and (m^k-1)/(m-1) are prime with\n";
   else fout << " is prime with\n";
   if (k > 1) fout << "k  = " << k << ", and\n";
   if (c1 < INT_MAX || c2 < INT_MAX) fout << "2^e + c1 < m < 2^e + c2, for e = " << e << ", c1 = "
         << c1 << ", c2 = " << c2 << "\n";
   else fout << "m < 2^e and e = " << e << "\n";
   fout << "Search parameters:\n";
   fout << "  safe = " << std::boolalpha << safe;
   fout << ",  factomm1 = " << factomm1;
   fout << ",  numtrials = " << numtrials << "\n\n";
}

} // end namespace LatMRG

#endif

// This file is part of LatMRG.
//
// Copyright (C) 2012-2023  The LatMRG authors, under the supervision
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
#include "latticetester/IntFactor.h"

#include "latmrg/Chrono.h"
#include "latmrg/IntFactorization.h"

namespace LatMRG {

  /**
   * This class provides methods to search for integers \f$m\f$ that are prime,
   * for which the integer \f$r = (m^k-1)/(m-1)\f$ is also prime for a given
   * \f$k\f$, and possibly for which \f$(m-1)/2\f$ is also prime.
   *
   * This class solely uses probabilistic primality tests, but does it a lot of
   * times.
   *
   * \todo Could probably switch to a deterministic test because it is possible
   * to do everything this class does with efficient deterministic tests for integers < 2^64
   * see https://math.stackexchange.com/questions/2481148/primality-testing-for-64-bit-numbers
   */
template<typename Int> static class PrimesFinder {

      public:

        /**
         * Constructor.
         */
        PrimesFinder();

        /**
         * Destructor.
         */
        ~PrimesFinder();

        /**
         * Finds \f$s\f$ prime integers \f$m<2^e\f$ that are closest to
         * \f$2^e\f$. If `facto` is `true`, then \f$m-1\f$ is factorized in its
         * prime factors. The results are printed on stream `fout`.
         */
        static void findPrime (int e, int s, bool facto, std::ostream & fout);

        /**
         * Finds \f$s\f$ prime integers \f$m<2^e\f$ that are closest to
         * \f$2^e\f$, such that \f$m\f$ and \f$r = (m^k-1)/(m-1)\f$ are prime.
         * If `safe` is `true`, \f$(m-1)/2\f$ is also required to be prime. The
         * results are printed on stream `fout`. If `facto` is `true`, then
         * \f$m-1\f$ is factorized in its prime factors. If \f$k=1\f$, \f$r\f$
         * is considered to be prime.
         */
        static void findPrime (int k, int e, int s, bool safe, bool facto, std::ostream & fout);

        /**
         * Finds all integers \f$m\f$, in \f$2^e + c_1 \le m \le2^e + c_2\f$,
         * such that \f$m\f$ and \f$r = (m^k-1)/(m-1)\f$ are prime. If `safe`
         * is `true`, \f$(m-1)/2\f$ is also required to be prime. The results
         * are printed on stream `fout`. If `facto` is `true`, then \f$m-1\f$
         * is factorized in prime factors. If \f$k=1\f$, \f$r\f$ is considered
         * to be prime.
         */
        static void findPrime (int k, int e, long c1, long c2, bool safe, bool facto,
            std::ostream & fout);

        /**
         * Used by the main executable to compare the examples.
         * */
        // bool print_time = true;

      private:

        /**
         * Writes the parameters of the find to the stream `fout`.
         */
        void writeHeader (int k, int e, long c1, long c2, bool safe, bool facto,
            std::ostream & fout);

        /**
         * Writes the CPU time of the find to the stream `fout`.
         */
        void writeFooter (std::ostream & fout);

        /**
         * This is the general purpose function used by this program. This
         * method searches for `s` prime integers between `S1` and `S2`.
         * */
        void findPrime (int k, int e, int s, const Int & S1, const Int & S2, bool safe,
            bool facto, std::ostream & fout, const long KTRIALS = 200);

        // Chrono timer;


        void nextM (Int & m) {
          m -= 2;
          if (0 == m % 5)
            m -= 2;
        }
    };

  //============================================================================

  template<typename Int>
    PrimesFinder<Int>::PrimesFinder ()
    {
      timer.init();
    }

  //===========================================================================

  template<typename Int>
    PrimesFinder<Int>::~PrimesFinder ()
    {
    }

  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::findPrime (int k, int e, int s, const Int & S1, const Int & S2,
        bool safe, bool facto, std::ostream & fout, const long KTRIALS)
    {
      Int m;
      if (NTL::IsOdd (S2))
        m = S2;
      else
        m = S2 - Int(1);
      int i = 0;

      while (i < s && m >= S1) {
        LatticeTester::PrimeType status = LatticeTester::IntFactor<Int>::isPrime (m, KTRIALS);
        if (status == LatticeTester::PRIME || status == LatticeTester::PROB_PRIME) {
          Int m1 = m - Int(1);
          if (safe) {
            if (1 == m % 4) {
              nextM(m);
              continue;
            }
            Int m1s2 = m1 / Int(2);
            status = LatticeTester::IntFactor<Int>::isPrime (m1s2, KTRIALS);
            if (status != LatticeTester::PRIME && status != LatticeTester::PROB_PRIME) {
              nextM(m);
              continue;
            }
          }
          Int r;
          NTL::set(r);
          if (k > 1) {
            r = NTL::power (m, k);
            --r;
            r = r/m1;
            status = LatticeTester::IntFactor<Int>::isPrime (r, KTRIALS);
          }
          if (k == 1 || status == LatticeTester::PRIME || status == LatticeTester::PROB_PRIME) {
            i++;
            fout << "//========================================================"
              "======================\n\n";
            fout << "m = " << m << "\n";
            Int Sdiff = m - (Int(1)<<e);
            fout << "  = 2^" << e;
            if (Sdiff >= 0) {
              fout << " + ";
              fout << Sdiff << "\n";
            } else {
              fout << " - ";
              fout << (-Sdiff) << "\n";
            }
            fout << "\n";
            if (facto) {
              IntFactorization<Int> ifac;
              ifac.clear();
              ifac.setNumber (m1);
              ifac.factorize ();
              fout << "Factors of m - 1 = ";
              fout << ifac.toString ();
            }
          }
        }
        nextM(m);
      }
    }

  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::findPrime (int k, int e, int s, bool safe, bool facto,
        std::ostream & fout)
    {
      Int Sm1, Sm2;
      writeHeader (k, e, INT_MAX, INT_MAX, safe, facto, fout);
      timer.init();
      Sm2 = (Int(1)<<e) - 1;
      Sm1 = 2;
      findPrime (k, e, s, Sm1, Sm2, safe, facto, fout);
      writeFooter (fout);
    }


  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::findPrime (int e, int s, bool facto, std::ostream & fout)
    {
      Int Sm1, Sm2;
      writeHeader (1, e, INT_MAX, INT_MAX, false, facto, fout);
      timer.init();
      Sm2 = (Int(1)<<e) - 1;
      Sm1 = 2;
      findPrime (1, e, s, Sm1, Sm2, false, facto, fout);
      writeFooter (fout);
    }


  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::findPrime (int k, int e, long c1, long c2, bool safe,
        bool facto, std::ostream & fout)
    {
      Int Sm1, Sm2;
      writeHeader (k, e, c1, c2, safe, facto, fout);
      timer.init();
      Sm1 = (Int(1)<<e) + c1;
      Sm2 = (Int(1)<<e) + c2;
      findPrime (k, e, INT_MAX, Sm1, Sm2, safe, facto, fout);
      writeFooter (fout);
    }


  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::writeHeader (int k, int e, long c1, long c2, bool safe,
        bool facto, std::ostream & fout)
    {
      fout << "-----------------------------------------------------" << std::endl;
      fout << "Values such that m";
      if (safe)
        fout << ", (m-1)/2,";
      if (k > 1)
        fout << " and (m^k-1)/(m-1) are prime with\n";
      else
        fout << " is prime with\n";
      fout << "k  = " << k << ", and\n";
      if (c1 < INT_MAX || c2 < INT_MAX)
        fout << "2^e + c1 < m < 2^e + c2, for e = " << e << " c1 = " << c1
          << " c2 = " << c2 << "\n";
      else
        fout << "m < 2^e, and e = " << e << "\n";

      fout << "\nProgram values:\n";
      fout << "safe = " << std::boolalpha << safe << "\n";
      fout << "facto = " << facto << "\n\n";
    }


  //===========================================================================

  template<typename Int>
    void PrimesFinder<Int>::writeFooter (std::ostream & fout)
    {
      fout << "\nCPU time: ";
      fout << timer.toString () << std::endl;
    }

}
#endif
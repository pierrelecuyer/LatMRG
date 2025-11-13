/**
 * TestMWCSpeed: compares speed of various MWC generators.  The order `k` goes from 1 to 3.
 * Some have a_0 = -1, others have a general a_0.  Some have forces zero coefficients a_j.
 */
#include <iostream>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <stdint.h>
#include <math.h>
#include <iomanip> // For std::hex, std::setfill, std::setw
#include <algorithm> // Required for std::reverse

uint64_t x, y, z, c;
uint64_t x1, x2, x3;
__uint128_t tau, sum3 = 0;
uint64_t sum = 0;
// State (x_{n-1},...,x_{n-k},c) = (x1,...,xk,c).
// State xx = (x_0,...,x_{k-1},c) = (xk,...,x1,c).

// Function to print a __uint128_t in decimal.
std::string uint128_dec_str(unsigned __int128 val) {
   if (val == 0) {
      return "0";
   }
   std::string s;
   while (val > 0) {
      s += (val % 10) + '0';
      val /= 10;
   }
   std::reverse(s.begin(), s.end());
   return s;
}

static void printState() {
   std::cout << "sum = " << sum << "\n";
   std::cout << "c   = " << c << "\n";
   std::cout << "tau   = " << uint128_dec_str(tau) << "\n";
}

static void printResults(const std::string &rngName, clock_t tmp, uint64_t sum) {
   std::cout << std::setw(16) << rngName << std::fixed << std::setw(13)
         << (double) tmp / (CLOCKS_PER_SEC) << "    " << std::setw(18) << sum << "\n";
}

static void printResultsDouble(const std::string &rngName, clock_t tmp, double average) {
   std::cout << std::setw(16) << rngName << std::setprecision(8) << std::fixed << std::setw(13)
         << (double) tmp / (CLOCKS_PER_SEC) << "  " << std::setw(12) << average << "\n";
}

// *******   k = 1  **********************************************

// From Vigna 2021
uint64_t inline MWC128() {
   const uint64_t result = x;
   const __uint128_t t = 0xffebb71d94fcdaf9 * (__uint128_t ) x + c;
   x = t;
   c = t >> 64;
   return result;
   // return x;
}

// k = 1, a_0 = -1.
uint64_t inline mwc64k1() {
   x2 = x1;
   tau = (0xffebb71d94fcdaf9 * (__uint128_t ) x1) + c;
   c = tau >> 64;
   x1 = tau;
   return x2;
}

// *******   k = 2  **********************************************

// From Vigna 2021.   Here, a_0=-1 and a_1=0.
uint64_t inline MWC192() {
   const uint64_t result = y;
   const __uint128_t tau = 0xffa04e67b3c95d86 * (__uint128_t ) x + c;
   x = y;
   y = tau;
   c = tau >> 64;
   return result;
}

// ****************************************************************
// From MWC1k2-0-62.res   Here, a_0=-1 and a_1=0.
// -0x1 0x0 0x7c309ee45ade
uint64_t inline mwc64k2one() {
   x3 = x1;
   tau = (0x7c309ee45ade * (__uint128_t ) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return x3;
}

// From MWC1k2-62-58.res   a_0=-1
// - 0x1 0x2ae390b92740f6d 0x6fcce264fcc37
uint64_t inline mwc64k2two() {
   const uint64_t out = x1;
   tau = (0x2ae390b92740f6d * (__uint128_t ) x1 + 0x6fcce264fcc37 * (__uint128_t ) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWCa0k2-0-62-58.res    a_0 general, a_1=0
// - 0xe4e1b5e445c2a65 0x0 0x6595a0b6737e
// inv = 0xba1485699989776d
uint64_t inline mwc64k2a0one() {
   tau = 0x6595a0b6737e * (__uint128_t ) x2 + c;
   x2 = x1;
   x1 = 0xba1485699989776d * (uint64_t) tau;
   c = (tau - 0xe4e1b5e445c2a65 * (__uint128_t ) x1) >> 64;
   return x2;
}

// From MWCa0k2-60-58.res   a_0 general
// - 0xb52048b94ffe40b 0x27612f107be7335 0x11716b6421c7e3
//  inv = 0xa206661c217087a3
uint64_t inline mwc64k2a0two() {
   const uint64_t out = x1;
   tau = (0x27612f107be7335 * (__uint128_t ) x1 + 0x11716b6421c7e3 * (__uint128_t ) x2) + c;
   x2 = x1;
   x1 = 0xa206661c217087a3 * (uint64_t) tau;
   c = (tau - 0xb52048b94ffe40b * (__uint128_t ) x1) >> 64;
   return out;
}

// *******   k = 3  *********************************************

// From Vigna 2021.    a_0=-1, a_1 = a_2 = 0.
uint64_t inline MWC256() {
   const uint64_t result = z;
   const __uint128_t tau = (0xfff62cf2ccc0cdaf * (__uint128_t ) x) + c;
   x = y;
   y = z;
   z = tau;
   c = tau >> 64;
   return result;
}

// Parameters from Vigna 2021.
#define GMWC_MINUSA0 0x54c3da46afb70f
#define GMWC_A0INV 0xbbf397e9a69da811
#define GMWC_A3 0xff963a86efd088a2

// From Vigna 2021.   a_0 general, a_1 = a_2 = 0.
uint64_t inline GMWC256() {
   const __uint128_t tau = GMWC_A3 * (__uint128_t ) x + c;
   x = y;
   y = z;
   z = GMWC_A0INV * (uint64_t) tau;
   c = (tau + GMWC_MINUSA0 * (__uint128_t ) z) >> 64;
   return z;
}

// ****************************************************************
// From MWC1k3-00-62.res    a_0 = -1, a_1 = a_2 = 0.
//  -0x1 0x0 0x0 0x19e24ee5997a
uint64_t inline mwc64k3one() {
   const uint64_t out = x1;
   tau = (0x19e24ee5997a * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   //return x1;
   return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
//  - 0x1 0x0 0x13294296084b09e 0x43dfe9959abe
uint64_t inline mwc64k3two() {
   const uint64_t out = x1;
   tau = (0x13294296084b09e * (__uint128_t ) x2 + 0x43dfe9959abe * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   //return x1;
   return out;
}

// From MWC1k3-60-58.res      a_0 = -1.
// - 0x1 0x2df69b79a87ec31 0x3c9789537bf2185 0xf8611733381c
uint64_t inline mwc64k3three() {
   const uint64_t out = x1;
   tau = (0x2df69b79a87ec31 * (__uint128_t ) x1 + 0x3c9789537bf2185 * (__uint128_t ) x2
         + 0xf8611733381c * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// Three coefficients are equal.  a_0 =-1,  a_1 = a_2 = a_3.
uint64_t inline mwc64k3equal() {
   const uint64_t out = x1;
   //const __uint128_t
   tau = 0x1C2A54C95D8F833B * sum3 + c;
   sum3 = (__uint128_t ) x1 + (__uint128_t ) x2 + (__uint128_t ) x3;
   // tau = (0x1C2A54C95D8F833B * ((__uint128_t ) x1 + (__uint128_t ) x2 + (__uint128_t ) x3)) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
uint64_t inline mwc64k3twoXor() {
   const uint64_t out = x1 ^ x3;
   tau = (0x13294296084b09e * (__uint128_t ) x2 + 0x43dfe9959abe * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-60-58.res      a_0 = -1.
// 0x25c95c0f8e357bf 0x3ffaf09909ce57a 0xaeced5e89bc93
uint64_t inline mwc64k3threeXor() {
   const uint64_t out = x1 ^ x2;
   tau = (0x2df69b79a87ec31 * (__uint128_t ) x1 + 0x3c9789537bf2185 * (__uint128_t ) x2
         + 0xf8611733381c * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// **************************************************************

// From MWCa0k3-00-62.res       a_0 general, a_1 = a_2 = 0.
// - 0xb0aa144f79aaced 0x0 0x0 0x2bc396ccf289
// inv = 0x3093058f991690e5
uint64_t inline mwc64k3a0one() {
   const uint64_t out = x1;
   // const __uint128_t
   tau = 0x2bc396ccf289 * (__uint128_t ) x3 + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x3093058f991690e5 * (uint64_t) tau;
   c = (tau - 0xb0aa144f79aaced * (__uint128_t ) x1) >> 64;
   return out;
}

// From MWCa0k3-0-60-58.res      a_0 general, a_1 = 0.
// - 0xec1c85e684758cd 0x0 0x32411d0d7063217 0x67235523e93e
// inv = 0xe39b835e23585405
uint64_t inline mwc64k3a0two() {
   const uint64_t out = x1;
   tau = (0x32411d0d7063217 * (__uint128_t ) x2 + 0x67235523e93e * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0xe39b835e23585405 * (uint64_t) tau;
   c = (tau - 0xec1c85e684758cd * (__uint128_t ) x1) >> 64;
   return out;
}

// From MWCa0k3-60-58.res      a_0 general.
// - 0xc9071286719f83b 0x453c969a880a60 0x3a79cceaaf5b155 0x114f867f357a7
// 0x9b6e12f406d620f3
uint64_t inline mwc64k3a0three() {
   const uint64_t out = x1;
   tau = (0x453c969a880a60 * (__uint128_t ) x1 + 0x3a79cceaaf5b155 * (__uint128_t ) x2
         + 0x114f867f357a7 * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x9b6e12f406d620f3 * (uint64_t) tau;
   c = (tau - 0xc9071286719f83b * (__uint128_t ) x1) >> 64;
   return out;
}

// *************************************************************************

/**
 * To test the speed of the generators, we could use the following, which takes the
 * RNG function as a parameter, but this slows down the loop considerably,
 * For this reason, we just repeat the code instead.
 */
typedef uint64_t (rngFunc)();

static void testLoop(const std::string &rngName, rngFunc func, int64_t n) {
   x = y = z = x1 = x2 = x3 = c = 12345;
   sum = 0;
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += func();
   }
   tmp = clock() - tmp;
   printResults(rngName, tmp, sum);
}

// *************************************************************************

int main() {
   // int64_t n = 4;
   // int64_t n = 1000*1000; // One million
   int64_t n = 1000 * 1000 * 1000; // One billion
   std::cout << "\n=============================================================\n";
   std::cout << "Time to generate n = " << n << " = " << std::scientific << (double) n
         << " numbers.\n";
   std::cout << "    Generator     Time (seconds)      Sum mod 2^{64} \n";
   clock_t tmp;      // To measure computing time.

   // *******   k = 1  *********************************************
   std::cout << "k = 1 \n";

   // Using this `testLoop` function is much simpler than the other code that follows
   // with all the repetitions, but it is much slower.
   testLoop("MWC128 as param", MWC128, n);
   testLoop("mwc64k1 as param", mwc64k1, n);
   std::cout << "\n";

   x = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC128();
   }
   tmp = clock() - tmp;
   printResults("MWC128", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k1();
   }
   tmp = clock() - tmp;
   printResults("mwc64k1", tmp, sum);

   // *******   k = 2  *********************************************
   std::cout << "k = 2 \n";

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC192();
   }
   tmp = clock() - tmp;
   printResults("MWC192", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2one();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2one", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2two();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2two", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a0one();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a0one", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a0two();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a0two", tmp, sum);

   // *******   k = 3  *********************************************
   std::cout << "k = 3 \n";

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC256();
   }
   tmp = clock() - tmp;
   printResults("MWC256", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3one();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3one", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3two();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3two", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3three();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3three", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3equal();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3equal", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3twoXor();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3twoXor", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3threeXor();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3threeXor", tmp, sum);
   std::cout << "\n";

   // ************************************************************
   // a_0 < -1

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += GMWC256();
   }
   tmp = clock() - tmp;
   printResults("GMWC256", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a0one();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a0one", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a0two();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a0two", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a0three();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a0three", tmp, sum);

   // ******************************************************
   std::cout << "\nUniform over (0,1) \n";

   x1 = x2 = x3 = c = 12345;
   const double twom53 = 1.0 / (double) ((uint64_t) 1 << 53);
   const double twom63 = 1.0 / (double) ((uint64_t) 1 << 63);
   const double twom64 = 0.5 / (double) ((uint64_t) 1 << 63);
   // onem53 = ldexp(1, -53);
   std::cout << "twom53 = " << std::hexfloat << std::setprecision(16) << twom53 << std::fixed
         << "\n";
   std::cout << "twom63 = " << std::hexfloat << std::setprecision(16) << twom63 << std::fixed
         << "\n";
   std::cout << "twom64 = " << std::hexfloat << std::setprecision(16) << twom64 << std::fixed
         << "\n";
   std::cout << "standard deviation = (4n)^{-1/2} =   " << std::setprecision(8)
         << 1.0 / sqrt(2.0 * n) << "\n";
   std::cout << "We repeat the tests three times.\n\n";
   std::cout << "    Generator        Time (seconds)   Average \n";
   double dsum = 0.0;
   x1 = x2 = x3 = c = 12345;

   for (int64_t r = 0; r < 3; r++) {

      dsum = 0.0;
      tmp = clock();
      for (int64_t i = 0; i < n; i++) {
         dsum += (mwc64k3two() >> 11) * twom53;
      }
      tmp = clock() - tmp;
      printResultsDouble("mwc64k3two U(0,1)   ", tmp, dsum / n);

      dsum = 0.0;
      tmp = clock();
      for (int64_t i = 0; i < n; i++) {
         dsum += (mwc64k3three() >> 11) * twom53;
      }
      tmp = clock() - tmp;
      printResultsDouble("mwc64k3three U(0,1) ", tmp, dsum / n);

      dsum = 0.0;
      tmp = clock();
      for (int64_t i = 0; i < n; i++) {
         dsum += (mwc64k3twoXor() >> 11) * twom53;
      }
      tmp = clock() - tmp;
      printResultsDouble("mwc64k3twoXor U(0,1)", tmp, dsum / n);

      dsum = 0.0;
      tmp = clock();
      for (int64_t i = 0; i < n; i++) {
         dsum += (mwc64k3two() >> 1) * twom63;
      }
      tmp = clock() - tmp;
      printResultsDouble("mwc64k3two U(0,1) 63", tmp, dsum / n);

      dsum = 0.0;
      tmp = clock();
      for (int64_t i = 0; i < n; i++) {
         dsum += (mwc64k3two()) * twom64;
      }
      tmp = clock() - tmp;
      printResultsDouble("mwc64k3two U(0,1) 64", tmp, dsum / n);
      std::cout << "\n";
   }

   return 0;
}

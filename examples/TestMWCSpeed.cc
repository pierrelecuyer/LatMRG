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
const double twom53 = 1.0 / (double) ((uint64_t) 1 << 53);
const double twom55 = 1.0 / (double) ((uint64_t) 1 << 55);
const double twom63 = 1.0 / (double) ((uint64_t) 1 << 63);
const double twom64 = 0.5 / (double) ((uint64_t) 1 << 63);
clock_t tmp;      // To measure computing time.
clock_t tottmp;   // Total computing time.

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

// Called at the beginning of a speed test.
static void inline init() {
   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
}

static void printState() {
   std::cout << "sum = " << sum << "\n";
   std::cout << "c   = " << c << "\n";
   std::cout << "tau   = " << uint128_dec_str(tau) << "\n";
}

static void inline printResults(const std::string &rngName, clock_t tmp, uint64_t sum) {
   // tmp = clock() - tmp;
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
   tau = (0x87eac24b8adc9 * (__uint128_t ) x1) + c;
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
uint64_t inline mwc64k2a1() {
   x3 = x1;
   tau = (0x4b740f53265d * (__uint128_t) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return x3;
}

// From MWC1k2-62-58.res   a_0=-1
uint64_t inline mwc64k2a2() {
   const uint64_t out = x1;
   tau = (0x2ae390b92740f6d * (__uint128_t) x1 + 0x6fcce264fcc37 * (__uint128_t) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k2-62-58.res   a_0=-1
uint64_t inline mwc64k2a2Xor() {
   const uint64_t out = x1 ^ x2;
   tau = (0x2ae390b92740f6d * (__uint128_t) x1 + 0x6fcce264fcc37 * (__uint128_t) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// a_0=-1
uint64_t inline mwc64k2a2eq() {
   const uint64_t out = x1;
   tau = 0xc45ec2462f86 * ((__uint128_t) x1 + (__uint128_t) x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// a_0=-1
uint64_t inline mwc64k2a2eq2() {
   const uint64_t out = x1;
   // Note: swapping these two lines would be wrong!
   sum3 = (__uint128_t) x1 + (__uint128_t) x2;
   tau = 0xc45ec2462f86 * sum3 + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWCa0k2-0-62-58.res    a_0 general, a_1=0
uint64_t inline mwc64k2a1gk() {
   tau = 0x6595a0b6737e * (__uint128_t ) x2 + c;
   x2 = x1;
   x1 = 0xba1485699989776d * (uint64_t) tau;
   c = (tau - 0xe4e1b5e445c2a65 * (__uint128_t ) x1) >> 64;
   return x2;
}

// From MWCa0k2-60-58.res   a_0 general
uint64_t inline mwc64k2a2gk() {
   const uint64_t out = x1;
   tau = (0x22ddf8545f51c3d * (__uint128_t) x1 + 0x3b4ba2cc0eb83d39 * (__uint128_t) x2) + c;
   x2 = x1;
   x1 = 0x5a9b2e257b3e1d95 * (uint64_t) tau;
   c = (tau - 0xf0cdf98a7bf45bd * (__uint128_t ) x1) >> 64;
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

// ****************************************************************
// From MWC1k3-00-62.res    a_0 = -1, a_1 = a_2 = 0.
uint64_t inline mwc64k3a1() {
   const uint64_t out = x1;
   tau = (0x14fabef33841d * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   //return x1;
   return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
uint64_t inline mwc64k3a2() {
   const uint64_t out = x1;
   tau = (0x2902ebc31ec3683 * (__uint128_t) x2 + 0x57baa090037 * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   //return x1;
   return out;
}

// From MWC1k3-60-58.res      a_0 = -1.
uint64_t inline mwc64k3a3() {
   const uint64_t out = x1;
   tau = (0x979f3660cbaca4 * (__uint128_t) x1 + 0x31dcb74b96510a8 * (__uint128_t) x2
   + 0x5194d4649dcd * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// Two coefficients are equal.  a_0 =-1,  a_2 = a_3.
uint64_t inline mwc64k3a2eq() {
   const uint64_t out = x1;
   sum3 = (__uint128_t ) x2 + (__uint128_t ) x3;
   tau = 0x1d0710107d5d * sum3 + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// Three coefficients are equal.  a_0 =-1,  a_1 = a_2 = a_3.
uint64_t inline mwc64k3a3eq() {
   const uint64_t out = x1;
   sum3 = (__uint128_t ) x1 + (__uint128_t ) x2 + (__uint128_t ) x3;
   tau = 0x1d495210185c * sum3 + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = a_2 = 0.
uint64_t inline mwc64k3a1Xor() {
   const uint64_t out = x1 ^ x3;
   tau = (0x14fabef33841d * (__uint128_t ) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
uint64_t inline mwc64k3a2Xor() {
   const uint64_t out = x1 ^ x3;
   tau = (0x2902ebc31ec3683 * (__uint128_t) x2 + 0x57baa090037 * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-60-58.res      a_0 = -1.
uint64_t inline mwc64k3a3Xor() {
   const uint64_t out = x1 ^ x2;
   tau = (0x979f3660cbaca4 * (__uint128_t) x1 + 0x31dcb74b96510a8 * (__uint128_t) x2
   + 0x5194d4649dcd * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// ****************************************************************
// This version never returns 0, uses a while
uint64_t inline mwc64k3a2No0w() {
   uint64_t out = mwc64k3a2();
   while (out == 0)
      out = mwc64k3a2();
   return out;
}

// This version never returns 0, uses a if.
uint64_t inline mwc64k3a2No0i() {
   uint64_t out = mwc64k3a2();
   if (out == 0) out = mwc64k3a2No0i();
   return out;
}

// This version never returns 0
uint64_t inline mwc64k3a2No0a() {
   uint64_t out = mwc64k3a2();
   out += (out == 0) * mwc64k3a2();
   return out;
}

// This version never returns 0, uses while.
uint64_t inline mwc64k3a3No0w() {
   uint64_t out = mwc64k3a3();
   while (out == 0)
      out = mwc64k3a3();
   return out;
}

// This version never returns 0.
uint64_t inline mwc64k3a3No0i() {
   uint64_t out = mwc64k3a3();
   if (out == 0) out = mwc64k3a3No0i();
   return out;
}

// This version never returns 0.
uint64_t inline mwc64k3a3No0a() {
   uint64_t out = mwc64k3a3();
   out += (out == 0) * mwc64k3a3();
   return out;
}

// **************************************************************

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

// From MWCa0k3-00-62.res       a_0 general, a_1 = a_2 = 0.
uint64_t inline mwc64k3a1gk() {
   const uint64_t out = x1;
   // const __uint128_t
   tau = 0xf53334560d * (__uint128_t ) x3 + c;
   x3 = x2;
   x2 = x1;
   x1 = 0xa55ad7c337410603 * (uint64_t) tau;
   c = (tau - 0x2ce8eef308854ab * (__uint128_t ) x1) >> 64;
   return out;
}

// From MWCa0k3-0-60-58.res      a_0 general, a_1 = 0.
uint64_t inline mwc64k3a2gk() {
   const uint64_t out = x1;
   tau = (0x1a43e76662f692c * (__uint128_t) x2 + 0x5a888e5764193 * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x9772ec11e3b050c5 * (uint64_t) tau;
   c = (tau - 0xe6f81c49aeeae0d * (__uint128_t ) x1) >> 64;
   return out;
}

// From MWCa0k3-60-58.res      a_0 general.
uint64_t inline mwc64k3a3gk() {
   const uint64_t out = x1;
   tau = (0x1427e7f2aedfa9c * (__uint128_t) x1 + 0x373ad6944cfb8d0 * (__uint128_t) x2
   + 0xf2125fd7e3997 * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0xdd3034dd040dab6b * (uint64_t) tau;
   c = (tau - 0xea08673f5e82943 * (__uint128_t ) x1) >> 64;
   return out;
}

// ****************************************************************
// U(0,1) generators

uint64_t block53;

// This version may return 0
double inline mwc64k2a2U01() {
   return (mwc64k2a2() >> 11) * twom53;
}

// This version never returns 0, uses a while
double inline mwc64k2a2U01w() {
   block53 = (mwc64k2a2() >> 11);
   while (block53 == 0)
      block53 = (mwc64k2a2() >> 11);
   return block53 * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k2a2U01i() {
   block53 = (mwc64k2a2() >> 11);
   if (block53 == 0) return mwc64k2a2U01i();
   return block53 * twom53;
}

double inline mwc64k2a2XorU01() {
   return (mwc64k2a2Xor() >> 11) * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k2a2XorU01i() {
   block53 = (mwc64k2a2Xor() >> 11);
   if (block53 == 0) return mwc64k2a2XorU01i();
   return block53 * twom53;
}

// This version may return 0
double inline mwc64k3a1U01() {
   return (mwc64k3a1() >> 11) * twom53;
}

// This version never returns 0, uses a while
double inline mwc64k3a1U01w() {
   block53 = (mwc64k3a1() >> 11);
   while (block53 == 0)
      block53 = (mwc64k3a1() >> 11);
   return block53 * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k3a1U01i() {
   block53 = (mwc64k3a1() >> 11);
   if (block53 == 0) return mwc64k3a1U01i();
   return block53 * twom53;
}

double inline mwc64k3a1XorU01() {
   return (mwc64k3a1Xor() >> 11) * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k3a1XorU01i() {
   block53 = (mwc64k3a1Xor() >> 11);
   if (block53 == 0) return mwc64k3a1XorU01i();
   return block53 * twom53;
}

// This version may return 0
double inline mwc64k3a2U01() {
   return (mwc64k3a2() >> 11) * twom53;
}

// This version never returns 0, uses a while
double inline mwc64k3a2U01w() {
   block53 = (mwc64k3a2() >> 11);
   while (block53 == 0)
      block53 = (mwc64k3a2() >> 11);
   return block53 * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k3a2U01i() {
   block53 = (mwc64k3a2() >> 11);
   if (block53 == 0) return mwc64k3a2U01i();
   return block53 * twom53;
}

double inline mwc64k3a2XorU01() {
   return (mwc64k3a2Xor() >> 11) * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k3a2XorU01i() {
   block53 = (mwc64k3a2Xor() >> 11);
   if (block53 == 0) return mwc64k3a2XorU01i();
   return block53 * twom53;
}

// This version never returns 0, uses a while
double inline mwc64k3a3U01w() {
   block53 = (mwc64k3a3() >> 11);
   while (block53 == 0)
      block53 = (mwc64k3a3() >> 11);
   return block53 * twom53;
}

// This version never returns 0, uses a if.
double inline mwc64k3a3U01i() {
   block53 = (mwc64k3a3() >> 11);
   if (block53 == 0) return mwc64k3a3U01i();
   return block53 * twom53;
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
   //int64_t n = 1000 * 1000 * (int64_t) 10000; // Ten billions
   std::cout << "\n=============================================================\n";
   std::cout << "Time to generate n = " << n << " = " << std::scientific << (double) n
         << " numbers.\n";
   std::cout << "    Generator     Time (seconds)      Sum mod 2^{64} \n";
   tottmp = clock();

   // *******   k = 1  *********************************************
   std::cout << "k = 1 \n";

   // Using this `testLoop` function is much simpler than the other code that follows
   // with all the repetitions, but it is much slower.
   testLoop("MWC128 given as a parameter to testLoop    ", MWC128, n);
   testLoop("mwc64k1 given as a parameter to testLoop   ", mwc64k1, n);
   testLoop("mwc64k3a2 given as a parameter to testLoop ", mwc64k3a2, n);
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
      sum += mwc64k2a1();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a1", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a2();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a2", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a2eq();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a2eq", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a2eq2();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a2eq2", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a1gk();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a1gk", tmp, sum);

   x1 = x2 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k2a2gk();
   }
   tmp = clock() - tmp;
   printResults("mwc64k2a2gk", tmp, sum);

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
      sum += mwc64k3a1();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a1", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a2();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a2", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a3();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a3", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a2eq();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a2eq", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a3eq();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a3eq", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a1Xor();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a1Xor", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a2Xor();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a2Xor", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a3Xor();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a3Xor", tmp, sum);
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
      sum += mwc64k3a1gk();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a1gk", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a2gk();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a2gk", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3a3gk();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3a3gk", tmp, sum);

   // ******************************************************
   std::cout << "\nUniform over (0,1) \n";
   // onem53 = ldexp(1, -53);
   /*
    std::cout << "twom53 = " << std::hexfloat << std::setprecision(16) << twom53 << std::fixed
    << "\n";
    std::cout << "twom63 = " << std::hexfloat << std::setprecision(16) << twom63 << std::fixed
    << "\n";
    std::cout << "twom64 = " << std::hexfloat << std::setprecision(16) << twom64 << std::fixed
    << "\n";
    */
   std::cout << "standard dev. for average = (4n)^{-1/2} =    " << std::setprecision(8)
         << 1.0 / sqrt(2.0 * n) << "\n";
   std::cout << "      Generator              Time (seconds)    Average \n";
   double dsum = 0.0;
   //   for (int64_t r = 0; r < 3; r++) {

   // k = 2

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k2a2U01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2U01    53           ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k2a2() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2 U(0,1) 63          ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k2a2() >> 11) * twom53 + twom55;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2 U(0,1) 53, + twom55", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k2a2U01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2U01i U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k2a2XorU01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2XorU01 U(0,1) 53    ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k2a2XorU01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2XorU01i U(0,1) 53   ", tmp, dsum / n);

   std::cout << "\n";

   // k = 3

   // Important: We do not want the block of 53 bits to be 0.
   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1 U(0,1) 53          ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a1U01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1U01 U(0,1) 53       ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a1U01w();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1U01w U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a1U01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1U01i U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a1XorU01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1XorU01 U(0,1) 53    ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a1XorU01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1XorU01i U(0,1) 53   ", tmp, dsum / n);

   // *************************************
   std::cout << "\n";

   // Important: We do not want the block of 53 bits to be 0.
   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2 U(0,1) 53          ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a2U01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2U01 U(0,1) 53       ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a2U01w();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2U01w U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a2U01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2U01i U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a2XorU01();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2XorU01 U(0,1) 53    ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a2XorU01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2XorU01i U(0,1) 53   ", tmp, dsum / n);

   // *************************************
   std::cout << "\n";

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a3() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3 U(0,1) 53          ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a3U01w();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3U01w U(0,1) 53      ", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a3U01i();
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3U01i U(0,1) 53      ", tmp, dsum / n);

   return 0;  // ******************************************************************

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k2a2()) * twom64;
   }
   // tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2 U(0,1) 64", clock() - tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k2a2Xor() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2Xor U(0,1) 53", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k2a2Xor() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2Xor U(0,1) 63", tmp, dsum / n);

   x1 = x2 = x3 = c = 12345;
   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k2a2Xor() * twom64;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k2a2Xor U(0,1) 64", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1 U(0,1) 53", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1 U(0,1) 63", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1()) * twom64;
   }
   // tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1 U(0,1) 64", clock() - tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1Xor() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1Xor U(0,1)", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1Xor() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1Xor U(0,1) 63", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a1Xor()) * twom64;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a1Xor U(0,1) 64", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2 U(0,1)   ", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2 U(0,1) 63", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2()) * twom64;
   }
   // tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2 U(0,1) 64", clock() - tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2Xor() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2Xor U(0,1) 53", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2Xor() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2Xor U(0,1) 63", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a2Xor()) * twom64;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a2Xor U(0,1) 64", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a3() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3 U(0,1) 53 ", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3a3() >> 1) * twom63;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3 U(0,1) 63 ", tmp, dsum / n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += mwc64k3a3() * twom64;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3a3 U(0,1) 64 ", tmp, dsum / n);

   std::cout << "\n";

   printResultsDouble("Total computing time: ", clock() - tottmp, 0);
   return 0;
}

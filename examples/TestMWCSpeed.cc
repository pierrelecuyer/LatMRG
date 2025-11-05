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
__uint128_t tau;
uint64_t sum = 0;

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

static void printResults(const std::string& rngName, clock_t tmp, uint64_t sum) {
   std::cout << std::setw(16) << rngName << std::fixed << std::setw(13) << (double) tmp / (CLOCKS_PER_SEC)
      << "    " << std::setw(18) << sum << "\n";
}

static void printResultsDouble(const std::string& rngName, clock_t tmp, double average) {
   std::cout << std::setw(16) << rngName << std::setprecision(8) << std::fixed << std::setw(13) << (double) tmp / (CLOCKS_PER_SEC)
      << "  " << std::setw(12) << average << "\n";
}



// *******   k = 1  **********************************************

// From Vigna 2021
uint64_t inline MWC128() {
   const uint64_t result = x;
   const __uint128_t t = 0xffebb71d94fcdaf9 * (__uint128_t)x + c;
   x = t;
   c = t >> 64;
   return result;
   // return x;
}

// k = 1, a_0 = -1.
uint64_t inline mwc64k1() {
   x2 = x1;
   tau = (0xffebb71d94fcdaf9 * (__uint128_t)x1) + c;
   c = tau >> 64;
   x1 = tau;
   return x2;
}

// *******   k = 2  **********************************************

// From Vigna 2021.   Here, a_0=-1 and a_1=0.
uint64_t inline MWC192() {
   const uint64_t result = y;
   const __uint128_t tau = 0xffa04e67b3c95d86 * (__uint128_t)x + c;
   x = y;
   y = tau;
   c = tau >> 64;
   return result;
}

// ****************************************************************
// From MWC1k2-0-62.res   Here, a_0=-1 and a_1=0.
// 0x1 0x0 0xca13a67c40562
uint64_t inline mwc64k2one() {
   x3 = x1;
   tau = (0xca13a67c40562 * (__uint128_t)x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return x3;
}

// From MWC1k2-62-60.res   a_0=-1
// 0x1 0x7b88c6ac008d039 0x1d4f74ad35355f
uint64_t inline mwc64k2two() {
   const uint64_t out = x1;
   // x3 = x1;
   // tau = (0x1d4f74ad35355f * (__uint128_t)x2) + c
   //      + 0x7b88c6ac008d039 * (__uint128_t)x1;
   tau = (0x7b88c6ac008d039 * (__uint128_t)x1
        + 0x1d4f74ad35355f * (__uint128_t)x2) + c;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWCa0k2-0-62-60.res    a_0 general, a_1=0
// 0xbd8d4633ccc5425 0x0 0xeb93d750b30b55
uint64_t inline mwc64k2a0one() {
   tau = 0xeb93d750b30b55 * (__uint128_t)x2 + c;
   x2 = x1;
   x1 = 0xcbbaaf9408bba7ad * (uint64_t)tau;
   c = (tau - 0xbd8d4633ccc5425 * (__uint128_t)x1) >> 64;
   return x2;
}

// From MWCa0k2-60-58.res   a_0 general
// 0xc4c1a707eaf2611 0x104f649a2f561bd 0x12664d383dcf24
uint64_t inline mwc64k2a0two() {
   const uint64_t out = x1;
   tau = (0x104f649a2f561bd * (__uint128_t)x1
        + 0x12664d383dcf24 * (__uint128_t)x2) + c;
   x2 = x1;
   x1 = 0xd88ad801b1188af1 * (uint64_t)tau;
   c = (tau - 0xc4c1a707eaf2611 * (__uint128_t)x1) >> 64;
   return out;
}


// *******   k = 3  *********************************************

// From Vigna 2021.    a_0=-1, a_1 = a_2 = 0.
uint64_t inline MWC256() {
   const uint64_t result = z;
   const __uint128_t tau = (0xfff62cf2ccc0cdaf * (__uint128_t)x) + c;
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
   const __uint128_t tau = GMWC_A3 * (__uint128_t)x + c;
   x = y;
   y = z;
   z = GMWC_A0INV * (uint64_t)tau;
   c = (tau + GMWC_MINUSA0 * (__uint128_t)z) >> 64;
   return z;
}

// ****************************************************************
// From MWC1k3-00-62.res    a_0 = -1, a_1 = a_2 = 0.
uint64_t inline mwc64k3one() {
   // const uint64_t out = x1;
   tau = (0xbe2654711cb2ef * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return x1;
   //return out;
}

// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
uint64_t inline mwc64k3two() {
   const uint64_t out = x1;
   tau = (0x320fbe97bef0f95 * (__uint128_t)x2
        + 0x4a1849ec18bfa6 * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// From MWC1k3-60-58.res      a_0 = -1.
// 0x25c95c0f8e357bf 0x3ffaf09909ce57a 0xaeced5e89bc93
uint64_t inline mwc64k3three() {
   const uint64_t out = x1;
   tau = (0x25c95c0f8e357bf  * (__uint128_t)x1
         + 0x3ffaf09909ce57a * (__uint128_t)x2
         + 0xaeced5e89bc93 * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// Three coefficients are equal.  a_0 =-1,  a_1 = a_2 = a_3.
uint64_t inline mwc64k3eq() {
   const uint64_t out = x1;
   //const __uint128_t
   tau = (0x1C2A54C95D8F833B * ((__uint128_t)x1 + (__uint128_t)x2 + (__uint128_t)x3)) + c;
   // const unsigned __int128 __uint128_t tau = a3 * (__uint128_t)x3 + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// **************************************************************

// From MWCa0k3-00-62.res       a_0 general, a_1 = a_2 = 0.
// 0xea8671bc56fd655 0x0 0x0 0x12a74a7bdb3e366
uint64_t inline mwc64k3a0one() {
   // const __uint128_t
   tau = 0x12a74a7bdb3e366 * (__uint128_t)x3 + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x5d521274e6f676fd * (uint64_t)tau;
   c = (tau - 0xea8671bc56fd655 * (__uint128_t)x1) >> 64;
   return x1;
}

// From MWCa0k3-0-60-58.res      a_0 general, a_1 = 0.
// 0xf120f2a68911d19 0x0 0xef13750cf589db 0x2330790338222e
uint64_t inline mwc64k3a0two() {
   tau = (0xef13750cf589db * (__uint128_t)x2
        + 0x2330790338222e * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x67e84aa4b658ef29 * (uint64_t)tau;
   c = (tau - 0xf120f2a68911d19 * (__uint128_t)x1) >> 64;
   return x1;
}

// From MWCa0k3-60-58.res      a_0 general.
// 0xbc558636b864a79 0x36af52ce713f04a 0x3d8ce87d76c3ce8 0xa872c5b1a6655
uint64_t inline mwc64k3a0three() {
   tau = (0x36af52ce713f04a * (__uint128_t)x1
        + 0x3d8ce87d76c3ce8 * (__uint128_t)x2
        + 0xa872c5b1a6655   * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = 0x3099e2fa7229ffc9 * (uint64_t)tau;
   c = (tau - 0xbc558636b864a79 * (__uint128_t)x1) >> 64;
   return x1;

}


// *************************************************************************

/**
 * To test the speed of the generators, we could use the following, which takes the
 * RNG function as a parameter, but this slows down the loop considerably,
 * For this reason, we just repeat the code instead.
 */
typedef uint64_t (rngFunc)();

static void testLoop(const std::string& rngName, rngFunc func, int64_t n) {
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
   int64_t n = 1000*1000*1000; // One billion
   std::cout << "\n=============================================================\n";
   std::cout << "Time to generate n = " << n << " = " << std::scientific << (double)n << " numbers.\n";
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
      sum += mwc64k3eq();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3eq", tmp, sum);

   // ************************************************************

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
   double twom53 = 1.0 / (double)((uint64_t)1 << 53);
   // onem53 = ldexp(1, -53);
   std::cout << "twom53 = " << std::hexfloat << std::setprecision(16) << twom53 << std::fixed << "\n";
   std::cout << "standard deviation = (4n)^{-1/2} =   " << std::setprecision(8) << 1.0 / sqrt(2.0 * n) << "\n\n";
   std::cout << "    Generator        Time (seconds)   Average \n";
   double dsum = 0.0;
   x1 = x2 = x3 = c = 12345;

   /*
   double u01;
   for (int64_t i = 0; i < 3; i++) {
      uint64_t rand64 = mwc64k3two();
      u01 = (rand64 >> 11) * twom53;
      dsum += u01;
      std::cout << "rand64 = " << rand64 << "\n";
      // std::cout << "1.0p-53 = " << 1.0p-53 << "\n";
      std::cout << std::fixed << "u01 = " << u01 << "\n";
   }
   std::cout << "\n";
   */

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3two() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3two U(0,1)   ", tmp, dsum/n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3two() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3two U(0,1)   ", tmp, dsum/n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3three() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3three U(0,1) ", tmp, dsum/n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3three() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3three U(0,1) ", tmp, dsum/n);

   dsum = 0.0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      dsum += (mwc64k3three() >> 11) * twom53;
   }
   tmp = clock() - tmp;
   printResultsDouble("mwc64k3three U(0,1) ", tmp, dsum/n);

   return 0;

}

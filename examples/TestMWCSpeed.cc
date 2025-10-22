/**
 * TestMWCSpeed: compares speed of various MWC generators.
 */
#include <iostream>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <iomanip> // For std::hex, std::setfill, std::setw
#include <algorithm> // Required for std::reverse

// Function to print __uint128_t in decimal (more complex)
std::string uint128_dec_str(unsigned __int128 val) {
    if (val == 0) {
        return "0";
    }
    std::string s;
    while (val > 0) {
        s += (val % 10) + '0';
        val /= 10;
    }
    std::reverse(s.begin(), s.end()); // Reverse to get correct order
    return s;
}

uint64_t x, y, z, c;
uint64_t x1, x2, x3;
__uint128_t tau;
uint64_t sum = 0;


#define MWC_A1 0xffebb71d94fcdaf9

uint64_t inline MWC128() {
   const uint64_t result = x;
   // const __uint128_t t = 0xffebb71d94fcdaf9 * (__uint128_t)x + c;
   const __uint128_t t = (MWC_A1 * (__uint128_t)x) + c;
   x = t;
   c = t >> 64;
   return result;
   // return x;
}

uint64_t inline mwc64k1() {
   const uint64_t out = x1;
   // const __uint128_t
   tau = (0xffebb71d94fcdaf9 * (__uint128_t)x1) + c;
   c = tau >> 64;
   x1 = tau;
   return out;
   // return x1;
}

#define MWC_A2 0xffa04e67b3c95d86

uint64_t inline MWC192() {
   const uint64_t result = y;
   const __uint128_t tau = MWC_A2 * (__uint128_t)x + c;
   x = y;
   y = tau;
   c = tau >> 64;
   return result;
}

uint64_t inline MWC256() {
   const uint64_t result = z;
   // const __uint128_tau
   tau = (0xfff62cf2ccc0cdaf * (__uint128_t)x) + c;
   x = y;
   y = z;
   z = tau;
   c = tau >> 64;
   return result;
}

uint64_t inline mwc64k3one() {
   // const uint64_t out = x1;
   // const __uint128_t
   // tau = (0xF33273FB0E428AB6 * (__uint128_t)x3) + c;
   tau = (0xfff62cf2ccc0cdaf * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return x1;
   //return out;
}

uint64_t inline mwc64k3two() {
   const uint64_t out = x1;
   //__uint128_t
   tau = (0x14d45f35679075f * (__uint128_t)x2
       + 0xcd7f660abc3db5bf * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

// 0xde53bdba779465 0x289479d488799a9 0xc7f9ce844f52dd6b
uint64_t inline mwc64k3three() {
   const uint64_t out = x1;
   //__uint128_t
   tau = (0xde53bdba779465  * (__uint128_t)x1
         + 0x289479d488799a9 * (__uint128_t)x2
         + 0xc7f9ce844f52dd6b * (__uint128_t)x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

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

#define GMWC_MINUSA0 0x54c3da46afb70f
#define GMWC_A0INV 0xbbf397e9a69da811
#define GMWC_A3 0xff963a86efd088a2

uint64_t inline GMWC256() {
   const __uint128_t tau = GMWC_A3 * (__uint128_t)x + c;
   x = y;
   y = z;
   z = GMWC_A0INV * (uint64_t)tau;
   c = (tau + GMWC_MINUSA0 * (__uint128_t)z) >> 64;
   return z;
}


static void printState() {
   std::cout << "sum = " << sum << "\n";
   std::cout << "c   = " << c << "\n";
   std::cout << "tau   = " << uint128_dec_str(tau) << "\n";
}

static void printResults(const std::string& rngName, clock_t tmp, uint64_t sum) {
   // std::cout << "\n=============================================================\n";
   std::cout << "Generator: " << rngName << "\n";
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n";
   std::cout << "Sum of these numbers, modulo 2^{64}: " << sum << "\n";
   std::cout << "=============================================================\n\n";
}

int main() {
   // int64_t n = 4;
   // int64_t n = 1000*1000; // One million
   int64_t n = 1000*1000*1000; // One billion
   std::cout << "\n=============================================================\n";
   std::cout << "Time to generate n = " << n << " numbers.\n";
   clock_t tmp;      // To measure computing time.

   x = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC128();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC128", tmp, sum);

   x = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC128();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC128", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k1();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k1", tmp, sum);

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC192();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC192", tmp, sum);

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC256();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC256", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3one();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3one", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3two();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3two", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3three();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3three", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3eq();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3eq", tmp, sum);

   x = y = z = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += GMWC256();
      //printState();
   }
   tmp = clock() - tmp;
   printResults("GMWC256", tmp, sum);

   return 0;

}

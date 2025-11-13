/**
 * TestMWCSpeed: compares speed of various MWC generators implemented in ZZ.
 * This is mostly used to check the correctness of uint128_t implementations.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"

Int b = NTL::power(Int(2), 64);
Int two128 = NTL::power(Int(2), 128);

Int x, y, z;
Int x1, x2, x3, c;
Int tau;
Int sum;

Int inline MWC128() {
   const Int result = x;
   tau = (NTL::conv<Int>(0xffebb71d94fcdaf9) * x) + c;
   x = tau % b;
   c = tau >> 64;
   return result;
}

Int inline mwc64k1() {
   const Int out = x1;
   tau = (NTL::conv<Int>(0xffebb71d94fcdaf9) * x1) + c;
   x1 = tau % b;
   c = tau >> 64;
   return out;
}

#define MWC_A2 NTL::conv<Int>(0xffa04e67b3c95d86)

Int inline MWC192() {
   const Int result = y;
   tau = (MWC_A2 * x) + c;
   x = y;
   y = tau % b;
   c = tau >> 64;
   return result;
}

const Int a3 = NTL::conv<Int>(0xfff62cf2ccc0cdaf);

Int inline MWC256() {
   const Int result = z;
   tau = (a3 * x) + c;
   x = y;
   y = z;
   z = tau % b;
   c = tau >> 64;
   return result;
}

Int inline mwc64k3one() {
   const Int out = x1;
   //const Int
   tau = (NTL::conv<Int>(0xfff62cf2ccc0cdaf) * x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau % b;
   c = tau >> 64;
   return out;
}

Int inline mwc64k3two() {
   const Int out = x1;
   tau = (NTL::conv<Int>(0x14d45f35679075f) * x2
        + NTL::conv<Int>(0xcd7f660abc3db5bf) * x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau % b;
   c = tau >> 64;
   return out;
}

Int inline mwc64k3three() {
   const Int out = x1;
   tau = (NTL::conv<Int>(0xde53bdba779465) * x1
        + NTL::conv<Int>(0x289479d488799a9) * x2
        + NTL::conv<Int>(0xc7f9ce844f52dd6b) * x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau % b;
   c = tau >> 64;
   return out;
}

Int inline mwc64k3eq() {
   const Int out = x1;
   tau = (NTL::conv<Int>(0x1C2A54C95D8F833B) * (x1 + x2 + x3)) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau % b;
   c = tau >> 64;
   return out;
}


static void printState() {
   std::cout << "sum = " << sum << "\n";
   std::cout << "c   = " << c << "\n";
   std::cout << "tau   = " << tau << "\n";
   if (tau > two128)
      std::cout << "***  tau exceeds 2^{128}  ***\n";
   if (sum > b)
      std::cout << "***  sum exceeds 2^{64}  ***\n";}

static void printResults(const std::string& rngName, clock_t tmp, Int sum) {
   // std::cout << "\n=============================================================\n";
   std::cout << "Generator: " << rngName << "\n";
   std::cout << "Total running time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n";
   std::cout << "Sum of these numbers, modulo 2^{64}: " << sum  << "\n";
   std::cout << "=============================================================\n";
}

int main() {
   //int64_t n = 4;
   //int64_t n = 1000*1000; // One million
   int64_t n = 1000*1000*1000; // One billion
   std::cout << "\n=============================================================\n";
   std::cout << "Time to generate n = " << n << " numbers.\n";
   clock_t tmp;      // To measure computing time.

   x = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC128();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC128", tmp, sum);

   x = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC128();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC128", tmp, sum);

   x1 = x2 = x3 = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k1();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k1", tmp, sum);

   x = y = z = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC192();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC192", tmp, sum);

   x = y = z = x1 = x2 = x3 = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += MWC256();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("MWC256", tmp, sum);

   x1 = x2 = x3 = c = (Int)12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3one();
      sum %= b;
      //printState();
   }
   tmp = clock() - tmp;
   printResults("mwc64k3one", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3two();
      sum %= b;
      //printState();
      if (tau > two128)
         std::cout << "***  tau exceeds 2^{128}  ***\n";
   }
   tmp = clock() - tmp;
   printResults("mwc64k3two", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3three();
      sum %= b;
      //printState();
      if (tau > two128)
         std::cout << "***  tau exceeds 2^{128}  ***\n";
   }
   tmp = clock() - tmp;
   printResults("mwc64k3three", tmp, sum);

   x1 = x2 = x3 = c = 12345;
   sum = 0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      sum += mwc64k3eq();
      sum %= b;
      //printState();
      if (tau > two128)
         std::cout << "***  tau exceeds 2^{128}  ***\n";
   }
   tmp = clock() - tmp;
   printResults("mwc64k3eq", tmp, sum);

   return 0;

}

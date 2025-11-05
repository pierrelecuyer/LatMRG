/**
 * TestMWCJump: tests speed for jump-ahead for various MWC generators.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <stdint.h>
#include <math.h>
#include <iomanip> // For std::hex, std::setfill, std::setw
#include <algorithm> // Required for std::reverse

#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"

using namespace LatMRG;

Int b = NTL::power(Int(2), 64);
Int m;
IntVec aa, xx;
Int a0inv, c;

static void printState(const std::string& rngName, IntVec &xx, Int &c) {
   std::cout << std::setw(16) << rngName << std::fixed
      << "  " << std::setw(4) << xx << std::setw(4) << c << "\n";
}

static void printResultsTime(const std::string& rngName, clock_t tmp, IntVec &xx, Int &c) {
   std::cout << std::setw(16) << rngName << std::fixed << std::setw(12) << (double) tmp / (CLOCKS_PER_SEC)
      << "  " << std::setw(4) << xx << std::setw(4) << c << "\n";
}

// *************************************************************************

/**
 * To test the correctness and the speed of the jump-ahead for MWC generators.
 */

static void testLoop(const std::string& rngName, Int &b, IntVec &aa, Int &a0inv, Int &jumpSize, int64_t n0, int64_t n) {
   clock_t tmp;      // To measure computing time.
   Int jumpMult, jumpMult2;
   int64_t k = aa.length()-1;
   xx.SetLength(k);
   for (int64_t j = 0; j < k; j++) xx[j] = 12345;
   c = 12345;
   m = computeLCGModulusMWC(b, aa);
   // jumpSize = conv<Int>(5000);
   getJumpMult(jumpMult, b, m, jumpSize);
   std::cout << "\n    Generator                       vector xx                       carry \n";
   std::cout << "   jumpSize = " << jumpSize << "\n";
   for (int64_t i = 0; i < n0; i++) {
      jumpAhead(xx, c, b, aa, m, a0inv, xx, c, jumpMult);
      printState(rngName, xx, c);
   }
   for (int64_t j = 0; j < k; j++) xx[j] = 12345;
   c = 12345;
   Int jumpSize2 = jumpSize * conv<Int>(n0);
   getJumpMult(jumpMult2, b, m, jumpSize2);
   jumpAhead(xx, c, b, aa, m, a0inv, xx, c, jumpMult2);
   std::cout << "   jumpSize2 = " << jumpSize2 << "\n";
   printState(rngName, xx, c);

   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      jumpAhead(xx, c, b, aa, m, a0inv, xx, c, jumpMult);
   }
   tmp = clock() - tmp;
   std::cout << "\n    Generator     Time (seconds)           vector xx                         carry \n";
   printResultsTime(rngName, tmp, xx, c);
}

// *************************************************************************

int main() {
   int64_t n0 = 4;
   int64_t n = 1000*1000; // One million
   //int64_t n = 1000*1000*1000; // One billion
   Int jumpSize = conv<Int>(5000);

   std::cout << "\n=============================================================\n";
   std::cout << "Test of correctness of jump-ahead for the MWC \n";
   std::cout << "  and time to make n = " << n << " steps ahead. \n";

   // *******   k = 2  *********************************************
   aa.SetLength(3);

   // ****************************************************************
   // mwc64k2two
   aa[0] = -1;
   aa[1] = 0x7b88c6ac008d039;
   aa[2] = 0x1d4f74ad35355f;
   a0inv = -1;
   testLoop("mwc64k2two", b, aa, a0inv, jumpSize, n0, n);

   // ****************************************************************
   // mwc64k2a0two
   aa[0] = -0xc4c1a707eaf2611;
   aa[1] = 0x104f649a2f561bd;
   aa[2] = 0x12664d383dcf24;
   a0inv = -0xd88ad801b1188af1;
   testLoop("mwc64k2a0two", b, aa, a0inv, jumpSize, n0, n);

   // *******   k = 3  *********************************************
   aa.SetLength(3);

   // ****************************************************************
   // mwc64k3three
   aa[0] = -1;
   aa[1] = 0x25c95c0f8e357bf;
   aa[2] = 0x3ffaf09909ce57a;
   aa[3] = 0xaeced5e89bc93;
   a0inv = -1;
   testLoop("mwc64k3three", b, aa, a0inv, jumpSize, n0, n);

   // ****************************************************************
   // mwc64k3a0three
   aa[0] = -0xbc558636b864a79;
   aa[1] = 0x36af52ce713f04a;
   aa[2] = 0x3d8ce87d76c3ce8;
   aa[3] = 0xa872c5b1a6655;
   a0inv = -0x3099e2fa7229ffc9;
   testLoop("mwc64k3a0three", b, aa, a0inv, jumpSize, n0, n);

   return 0;

}

/**
 * TestMWCJump: tests speed for jump-ahead for various MWC generators.
 * The jump ahead is made with the LCG representation, using the NTL::ZZ type.
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

#include <gmp.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"

using namespace LatMRG;

Int b = NTL::power(Int(2), 64);
Int m;
IntVec aa, xx;
Int a0inv;
Int yy;
Int jumpMultHard;
Int jumpSizeHard = conv<Int>(1) << 60;  // 2^{60}
uint64_t x0, x1, x2, x3, c;
mpz_t bmp, mmp, ymp, jumpMultmp, sigmamp, xx0, xx1, xx2, xx3;


static void printState(IntVec &xx) {
   std::cout << std::fixed << "  " << std::setw(4) << xx << "\n";
}

static void printResultsTime(IntVec &xx, clock_t tmp) {
   std::cout << "  " << std::setw(4) << xx << std::fixed << std::setw(10)
        << (double) tmp / (CLOCKS_PER_SEC) << "\n";
}

// *************************************************************************

/**
 * To test the correctness and the speed of the jump-ahead for MWC generators.
 */
static void testLoop(const std::string& rngName, Int &jumpSize, int64_t n0, int64_t n) {
   clock_t tmp;      // To measure computing time.
   Int jumpMult, jumpMult2;
   int64_t k = aa.length()-1;
   xx.SetLength(k+1);
   std::cout << "\n=======================================================================\n";
   std::cout << "Generator name: " << rngName << ", k = " << k << ", a_0 = " << aa[0] << "\n";

   for (int64_t j = 0; j <= k; j++) xx[j] = 12345;
   m = computeLCGModulusMWC(b, aa);
   // jumpSize = conv<Int>(5000);
   getJumpMult(jumpMult, b, m, jumpSize);
   std::cout << "\n                               vector xx \n";
   std::cout << "  jumpSize = " << jumpSize << ",  numb jumps = " << n0 << "\n";
   for (int64_t i = 0; i < n0; i++) {
      jumpAhead(xx, b, aa, m, a0inv, xx, jumpMult);
      printState(xx);
   }
   for (int64_t j = 0; j < k; j++) xx[j] = 12345;
   Int jumpSize2 = jumpSize * conv<Int>(n0);
   getJumpMult(jumpMult2, b, m, jumpSize2);
   jumpAhead(xx, b, aa, m, a0inv, xx, jumpMult2);
   std::cout << "  jumpSize2 = " << jumpSize2 << ",  numb jumps = 1 \n";
   printState(xx);

   std::cout << "\n  Time to make n = " << n << " jumps ahead, from xx to next xx, "
         "jump size = " << jumpSize << "\n";
   std::cout << "                   last vector xx                            time (seconds) \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = 12345;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      jumpAhead(xx, b, aa, m, a0inv, xx, jumpMult);
   }
   tmp = clock() - tmp;
   printResultsTime(xx, tmp);

   std::cout << "\n  Jumps from LCG state y + transform to xx state: \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = 12345;
   Int y;
   MWCtoLCGState(y, b, aa, xx);
   std::cout << "  state y = " << y << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(y, jumpMult, y, m);
      LCGtoMWCState(xx, b, aa, a0inv, y);
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << y << "\n";
   printResultsTime(xx, tmp);
}


#define a1 0
#define a2 0x320fbe97bef0f95
#define a3 0x4a1849ec18bfa6

// Transform LCG state `yy` to MWC state `xx`, for k=3 and a_0=-1.
void inline LCGtoMWCHardk3Three1() {
   Int sigma = yy;
   xx[0] = sigma  & (b-1);
   // x3 = conv<ulong>(sigma);
   sigma = (sigma - xx[0]) >> 64;
   xx[1] = (sigma + a1 * xx[0]) & (b-1);
   sigma = (sigma - xx[1] + a1 * xx[0]) >> 64;
   xx[2] = (sigma + a1 * xx[1] + a2 * xx[0]) & (b-1);
   sigma = (sigma - x2 + a1 * xx[1] + a2 * xx[0]) >> 64;
   xx[3] = sigma;
}

// Transform LCG state `yy` to MWC state `xx`, for k=3, a_0=-1 and a_1=0.
void inline LCGtoMWCHardk3Two1() {
   Int sigma = yy;
   xx[0] = sigma  & (b-1);
   // x3 = conv<ulong>(sigma);
   sigma = (sigma - xx[0]) >> 64;
   xx[1] = (sigma) & (b-1);
   sigma = (sigma - xx[1]) >> 64;
   xx[2] = (sigma + a2 * xx[0]) & (b-1);
   xx[3] = (sigma - x2 + a2 * xx[0]) >> 64;
   x0 = conv<ulong>(xx[0]);
   x1 = conv<ulong>(xx[1]);
   x2 = conv<ulong>(xx[2]);
   x3 = conv<ulong>(xx[3]);
}

// Transform LCG state `yy` to MWC state `x0,...,x3`, for k=3, a_0=-1 and a_1=0.
void inline LCGtoMWCHardk3Two3() {
   Int sigma = yy;
   x0 = conv<ulong>(sigma);
   sigma = (sigma - conv<Int>(x0)) >> 64;
   x1 = conv<ulong>(sigma);
   sigma = (sigma - conv<Int>(x1)) >> 64;
   x2 = conv<ulong>(sigma + a2 * conv<Int>(x0));
   sigma = (sigma - conv<Int>(x2) + a2 * conv<Int>(x0)) >> 64;
   x3 = conv<ulong>(sigma);
}

/*
// Transform LCG state `ymp` to MWC state `x0,...,x3`, for k=3, a_0=-1 and a_1=0.
void inline LCGtoMWCHardk3Two4() {
   sigmamp = ymp;
   x0 = conv<ulong>(sigmamp);
   sigmamp = (sigmamp - conv<mpz_t>(x0)) >> 64;
   x1 = conv<ulong>(sigmamp);
   sigmamp = (sigmamp - conv<mpz_t>(x1)) >> 64;
   x2 = conv<ulong>(sigmamp + a2 * conv<mpz_t>(x0));
   sigmamp = (sigmamp - conv<mpz_t>(x2) + a2 * conv<mpz_t>(x0)) >> 64;
   x3 = conv<ulong>(sigmamp);
}
*/

// *************************************************************************

int main() {
   int64_t n0 = 4;
   int64_t n = 1000*1000; // One million
   //int64_t n = 1000*1000*1000; // One billion
   Int jumpSize = conv<Int>(5000);
   clock_t tmp;      // To measure computing time.

   std::cout << "\n=============================================================\n";
   std::cout << "Test of correctness of jump-ahead speed for the MWC \n";
   std::cout << "  and time to make n = " << n << " steps ahead. \n";
   std::cout << "The jump ahead is made using the LCG representation, using NTL::ZZ.\n";

   // *******   k = 2  *********************************************
   aa.SetLength(3);

   // ****************************************************************
   // mwc64k2two
   aa[0] = -1;
   aa[1] = 0x7b88c6ac008d039;
   aa[2] = 0x1d4f74ad35355f;
   a0inv = -1;
   testLoop("mwc64k2two", jumpSize, n0, n);

   // ****************************************************************
   // mwc64k2a0two
   aa[0] = -0xc4c1a707eaf2611;
   aa[1] = 0x104f649a2f561bd;
   aa[2] = 0x12664d383dcf24;
   a0inv = -0xd88ad801b1188af1;
   testLoop("mwc64k2a0two", jumpSize, n0, n);

   // *******   k = 3  *********************************************
   aa.SetLength(4);

   // ****************************************************************
   // mwc64k3two
   aa[0] = -1;
   aa[1] = 0;
   aa[2] = 0x320fbe97bef0f95;
   aa[3] = 0x4a1849ec18bfa6;
   a0inv = -1;
   testLoop("mwc64k3two", jumpSize, n0, n);

   /*
   // ****************************************************************
   // mwc64k3three
   aa[0] = -1;
   aa[1] = 0x25c95c0f8e357bf;
   aa[2] = 0x3ffaf09909ce57a;
   aa[3] = 0xaeced5e89bc93;
   a0inv = -1;
   testLoop("mwc64k3three", jumpSize, n0, n);

   // ****************************************************************
   // mwc64k3a0three
   aa[0] = -0xbc558636b864a79;
   aa[1] = 0x36af52ce713f04a;
   aa[2] = 0x3d8ce87d76c3ce8;
   aa[3] = 0xa872c5b1a6655;
   a0inv = -0x3099e2fa7229ffc9;
   testLoop("mwc64k3a0three", jumpSize, n0, n);
*/

   // ****************************************************************
   // mwc64k3two
   std::cout << "\n=======================================================================\n";
   std::cout << "\n Generator name: mcw64k3two, hardcoded jumps from LCG state y only, no transform to xx. \n";
   std::cout << "\n  Time to make n = " << n << " jumps ahead, jump size = " << jumpSize << "\n";
   getJumpMult(jumpMultHard, b, m, jumpSize);
   // x1 = x2 = x3 = c = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGState(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << yy << "\n";
   LCGtoMWCHardk3Three1();
   printResultsTime(xx, tmp);
   //std::cout << "  " << std::setw(4) << x3 << " "  << x2 << " " << x1 << " "  << c << " "
   //     << std::fixed << std::setw(10) << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";

   std::cout << "\n Generator name: mcw64k3two, hardcoded jumps from LCG state y + toMWC. \n";
   std::cout << "\n  Time to make n = " << n << " jumps ahead, jump size = " << jumpSize << "\n";
   // getJumpMult(jumpMultHard, b, m, jumpSize);
   // x1 = x2 = x3 = c = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGState(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3Three1();
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << yy << "\n";
   printResultsTime(xx, tmp);
   //std::cout << "  " << std::setw(4) << x3 << " "  << x2 << " " << x1 << " "  << c << " "
   //     << std::fixed << std::setw(10) << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";

   std::cout << "\n Generator name: mcw64k3two, hardcoded jumps from LCG state y, force a_1 = 0. \n";
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGState(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3Two1();
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << yy << "\n";
   printResultsTime(xx, tmp);

   std::cout << "\n Generator name: mcw64k3two, hardcoded jumps from LCG state y + toMWC. \n";
   std::cout << "\n  Time to make n = " << n << " jumps ahead, jump size = " << jumpSize << "\n";
   // getJumpMult(jumpMultHard, b, m, jumpSize);
   x0 = x1 = x2 = x3 = c = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGState(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3Two3();
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << yy << "\n";
   // printResultsTime(xx, tmp);
   std::cout << "      " << std::setw(4) << x0 << " "  << x1 << " " << x2 << " "  << x3 << " "
        << std::fixed << std::setw(10) << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";

   /*
   // mpz_t bmp, mmp, ymp, jumpMultmp, sigmamp, xx0, xx1, xx2, xx3;
   std::cout << "\n Generator name: mcw64k3two, jumps from LCG state y + toMWC, using GMP. \n";
   std::cout << "\n  Time to make n = " << n << " jumps ahead, jump size = " << jumpSize << "\n";
   // getJumpMult(jumpMultHard, b, m, jumpSize);
   jumpMultmp = conv<mpz_t>(jumpMultHard);
   xx0 = xx1 = xx2 = xx3 = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGState(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   ymp = conv<mpz_t>(yy);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      ymp = conv<mpz_t>(yy);
      LCGtoMWCHardk3Two4();
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << ymp << "\n";
   // printResultsTime(xx, tmp);
   std::cout << "      " << std::setw(4) << xx0 << " "  << xx1 << " " << xx2 << " "  << xx3 << " "
        << std::fixed << std::setw(10) << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
  */

   return 0;

}

/**
 * TestMWCJump: correctness and speed tests for jump-ahead for various MWC generators.
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
// #include "latmrg/StreamMWC.h"

using namespace LatMRG;

Int b = NTL::power(Int(2), 64);
Int m;
IntVec aa, xx;  //
Int a0inv;
Int yy, yy2;
Int jumpMultHard;
Int jumpSizeHard = conv<Int>(1) << 60;  // 2^{60}
uint64_t x0, x1, x2, x3, c;  // Here we assume that these integers cannot be negative!
// mpz_t bmp, mmp, ymp, jumpMultmp, sigmamp, xx0, xx1, xx2, xx3;


static void printStatexx() {
   std::cout << "   MWC state = " << std::setw(4) << xx << std::fixed << std::setw(10) << "\n";
}

static void printStatex() {
   std::cout << "   MWC state = " << std::setw(4) << x0 << " "  << x1 << " " << x2 << " "  << x3 << "\n";
}

static void printResultsTimex(int64_t n, Int &jumpSize, clock_t tmp) {
   std::cout << "   Time to make n = " << n << " jumps ahead, jump size = " << jumpSize
          << ",     CPU time:  " << (double) tmp / (CLOCKS_PER_SEC) << "\n";
   std::cout << "   final LCG state y = " << yy << "\n";
   printStatex();
}

static void printResultsTimexx(int64_t n, Int &jumpSize, clock_t tmp) {
   std::cout << "   Time to make n = " << n << " jumps ahead, jump size = " << jumpSize
          << ",     CPU time:  " << (double) tmp / (CLOCKS_PER_SEC) << "\n";
   std::cout << "   final LCG state y = " << yy << "\n";
   printStatexx();
}

/**
 * Prints `nsucc` successive states when starting from LCG state y=1,
 * first when running backward, then when running forward.
 */
static void succStates(int64_t nsucc) {
   std::cout << "  Successive states when starting from LCG state y=1, \n";
   std::cout << "  first when running backward, then when running forward. \n";
   std::cout << "    state-y   xx       y   y_{n-k}    xx   next-y \n";
   Int mult = b;  // The LCG multiplier.
   for (int64_t dir = 0; dir < 2; dir++) {
      if (dir == 1) mult = computeInvb(m, b);
      std::cout << "    mult = " << mult << "\n";
      yy = 1;
      for (int64_t i = 0; i < nsucc; i++) {
         LCGtoMWCStateDirect(xx, b, aa, m, yy);
         std::cout << "    " << std::setw(4) << yy << "     " << xx << "     ";
         MWCtoLCGState(yy2, b, aa, xx);
         std::cout << yy2 << "     ";
         MWCtoLCGStateLagk(yy2, b, aa, xx);
         std::cout << yy2 << "     ";
         LCGtoMWCStateLagk(xx, b, aa, a0inv, yy2);
         std::cout << xx << "     ";
         NTL::MulMod(yy, mult, yy, m);  // Move LCG by one step.
         std::cout << yy << "\n";
      }
   std::cout << "\n";
   }
}

 /**
 * Tests the correctness and speed of the jump-ahead for an MWC generator.
 * We perform `n0` jumps of size `jumpSize`, then one jump of size `n0 * jumpSize`
 * from the same initial state for comparison, then `n` jumps of size `jumpSize`
 * for the MWC state for timing, then `n` jumps of size `jumpSize` with the LCG
 * state with a transformation to the MWC state after each jump, also for timing.
 */
static void testLoop(const std::string& rngName, Int &jumpSize, int64_t nsucc, int64_t n0, int64_t n, int64_t state0) {
   clock_t tmp;      // To measure computing time.
   Int jumpMult, jumpMult2;  // Multipliers for jumping ahead with the LCG.
   int64_t k = aa.length()-1;
   xx.SetLength(k+1);        // Vector to save the MWC state.
   std::cout << "\n=======================================================================\n";
   std::cout << "Generator name: " << rngName << ", k = " << k << ", a_0 = " << aa[0] << "\n";
   m = computeLCGModulusMWC(b, aa);
   std::cout << "  m = " << m << "\n\n";

   if (nsucc > 0) succStates(nsucc);

   // Perform `n0` jumps of size `jumpSize`.
   getJumpAheadMult(jumpMult, b, m, jumpSize);
   std::cout << "  Perform successive jumps ahead and print MWC state after each jump.\n";
   std::cout << "  jumpSize = " << jumpSize << ",  numb jumps = " << n0 << "\n";
   std::cout << "  With jumpMWC:\n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;  // Same initial state for all experiments.
   for (int64_t i = 0; i < n0; i++) {
      jumpMWC(xx, b, aa, m, a0inv, xx, jumpMult, yy);
      printStatexx();
      // std::cout << "   state y_{n-k} = " << yy << "\n";
   }

   // Perform `n0` jumps of size `jumpSize` via the direct method.
   std::cout << "  With jumpMWCDirect:\n";
   // std::cout << "  jumpSize = " << jumpSize << ",  numb jumps = " << n0 << "\n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;  // Same initial state for all experiments.
   for (int64_t i = 0; i < n0; i++) {
      jumpMWCDirect(xx, b, aa, m, xx, jumpMult, yy);
      printStatexx();
      // std::cout << "   state y_{n} = " << yy << "\n";
   }

   // Perform one jump of size `n0 * jumpSize`
   std::cout << "\n  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize \n";
   Int jumpSize2 = jumpSize * conv<Int>(n0);
   getJumpAheadMult(jumpMult2, b, m, jumpSize2);
   std::cout << "  jumpSize2 = " << jumpSize2 << ",  numb jumps = 1, MWC state after the jump: \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;
   jumpMWC(xx, b, aa, m, a0inv, xx, jumpMult2, yy);
   printStatexx();

   // The following are timing tests. They also test if the state is correct after the jumps.
   // Perform `n` jumps of size `jumpSize` for the MWC via the LCG, using `jumpMWC`.
   std::cout << "\n  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      jumpMWC(xx, b, aa, m, a0inv, xx, jumpMult, yy);
   }
   tmp = clock() - tmp;
   printResultsTimexx(n, jumpSize, tmp);

   // Perform `n` jumps of size `jumpSize` for the LCG and compute the MWC state after each jump.
   std::cout << "\n  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMult, yy, m);
      LCGtoMWCStateLagk(xx, b, aa, a0inv, yy);
   }
   tmp = clock() - tmp;
   printResultsTimexx(n, jumpSize, tmp);

   // Perform `n` jumps of size `jumpSize` for the LCG and compute the MWC state after each jump,
   // this time by using LCGtoMWCStateDirect.
   std::cout << "\n  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. \n";
   for (int64_t j = 0; j <= k; j++) xx[j] = state0;
   MWCtoLCGState(yy, b, aa, xx);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMult, yy, m);
      LCGtoMWCStateDirect(xx, b, aa, m, yy);
   }
   tmp = clock() - tmp;
   printResultsTimexx(n, jumpSize, tmp);
}


/**
 * In the following, we define specialized functions to transform the LCG state
 * to a MWC state, for k=3 and 2 or 3 nonzero coefficients a_j.
 * We want to compare their performance with the general implementation.
 */
#define a1 0
#define a2 0x320fbe97bef0f95
#define a3 0x4a1849ec18bfa6

// Transform LCG state `yy` to MWC state `xx`, for k=3 and a_0=-1.
void inline LCGtoMWCHardk3a3() {
   Int sigma = yy;
   xx[0] = sigma  & (b-1);
   sigma = (sigma - xx[0]) >> 64;
   xx[1] = (sigma + a1 * xx[0]) & (b-1);
   sigma = (sigma - xx[1] + a1 * xx[0]) >> 64;
   xx[2] = (sigma + a1 * xx[1] + a2 * xx[0]) & (b-1);
   sigma = (sigma - x2 + a1 * xx[1] + a2 * xx[0]) >> 64;
   xx[3] = sigma;
}

// Transform LCG state `yy` to MWC state `xx`, for k=3, a_0=-1 and a_1=0.
// After that, `xx` is transformed to `uint64_t` variables `x0,...,x3`.
// This code assumes that `c_n` is never negative!!!
void inline LCGtoMWCHardk3a2xx() {
   Int sigma = yy;
   xx[0] = sigma  & (b-1);
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
void inline LCGtoMWCHardk3a2() {
   Int sigma = yy;
   x0 = conv<ulong>(sigma);
   sigma = (sigma - conv<Int>(x0)) >> 64;
   x1 = conv<ulong>(sigma);
   sigma = (sigma - conv<Int>(x1)) >> 64;
   x2 = conv<ulong>(sigma + a2 * conv<Int>(x0));
   sigma = (sigma - conv<Int>(x2) + a2 * conv<Int>(x0)) >> 64;
   x3 = conv<ulong>(sigma);
}

// Transform LCG state `yy` to MWC state `x0,...,x3`, for k=3, a_0=-1 and a_1=0.
// Here, some useless subtractions are removed.
void inline LCGtoMWCHardk3a2s() {
   Int sigma = yy;
   x0 = conv<ulong>(sigma);
   sigma = sigma >> 64;
   x1 = conv<ulong>(sigma);
   sigma = (sigma >> 64) + a2 * conv<Int>(x0);
   x2 = conv<ulong>(sigma);
   sigma = sigma >> 64;
   x3 = conv<ulong>(sigma);
}

/*
// Transform LCG state `ymp` to MWC state `x0,...,x3`, for k=3, a_0=-1 and a_1=0, using GMP.
void inline LCGtoMWCHardk3a2c() {
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
   int64_t n0 = 5;
   //int64_t n = 4; // One million
   int64_t n = 1000*1000; // One million
   //int64_t n = 1000*1000*1000; // One billion
   Int jumpSize = conv<Int>(5000);
   clock_t tmp;      // To measure computing time.

   std::cout << "\n=============================================================\n";
   std::cout << "Test of correctness and speed of jump-ahead or the MWC. \n";
   std::cout << "The jump ahead is made with the LCG representation, using NTL::ZZ.\n";


   // *******   k = 1  *********************************************

   // ****************************************************************
   // Baby case with k=1.
   aa.SetLength(2);
   aa[0] = -1;
   aa[1] = 3;
   a0inv = -1;
   b = NTL::power(Int(2), 2);
   jumpSize = conv<Int>(3);
   testLoop("baby11, k=1, b=4, m=11", jumpSize, 8, n0, n, 1);


   // *******   k = 2  *********************************************
   aa.SetLength(3);

   aa[0] = -1;
   aa[1] = 0;
   aa[2] = 3;
   a0inv = -1;
   // b = NTL::power(Int(2), 2);
   // jumpSize = conv<Int>(3);
   testLoop("baby47, k=2, b=4, m=47", jumpSize, 16, n0, n, 1);

   // return 0;

   // ****************************************************************
   b = NTL::power(Int(2), 64);
   jumpSize = conv<Int>(5000);
   n0 = 5;

   // ****************************************************************
   // mwc64k2a2
   aa[0] = -1;
   aa[1] = 0x7b88c6ac008d039;
   aa[2] = 0x1d4f74ad35355f;
   a0inv = -1;
   testLoop("mwc64k2a2", jumpSize, 0, n0, n, 12345);

   // return 0;

   // ****************************************************************
   // mwc64k2a2gk
   aa[0] = -0xc4c1a707eaf2611;
   aa[1] = 0x104f649a2f561bd;
   aa[2] = 0x12664d383dcf24;
   a0inv = -0xd88ad801b1188af1;
   testLoop("mwc64k2a2gk", jumpSize, 0, n0, n, 12345);


   // *******   k = 3  *********************************************
   aa.SetLength(4);

   // ****************************************************************
   // mwc64k3a2
   aa[0] = -1;
   aa[1] = 0;
   aa[2] = 0x320fbe97bef0f95;
   aa[3] = 0x4a1849ec18bfa6;
   a0inv = -1;
   testLoop("mwc64k3a2", jumpSize, 0, n0, n, 12345);

   // return 0;

   /*
   // ****************************************************************
   // mwc64k3three
   aa[0] = -1;
   aa[1] = 0x25c95c0f8e357bf;
   aa[2] = 0x3ffaf09909ce57a;
   aa[3] = 0xaeced5e89bc93;
   a0inv = -1;
   testLoop("mwc64k3three", jumpSize, n0, n);
   */

   // ****************************************************************
   // Here we test the specialized LCGtoMWC implementations to compare the speeds.

   // mwc64k3a2
   std::cout << "\n=======================================================================\n";
   std::cout << "\nGenerator name: mcw64k3a3, hardcoded jumps from LCG state y only, no transform to xx. \n";
   getJumpAheadMult(jumpMultHard, b, m, jumpSize);
   // x1 = x2 = x3 = c = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   // std::cout << "  state y = " << yy << "\n";
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
   }
   tmp = clock() - tmp;
   LCGtoMWCHardk3a3();
   printResultsTimexx(n, jumpSize, tmp);

   std::cout << "\nGenerator name: mcw64k3a2xx, hardcoded jumps from LCG state y + toMWC, using xx. \n";
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3a2xx();
   }
   tmp = clock() - tmp;
   printResultsTimexx(n, jumpSize, tmp);

   std::cout << "\nGenerator name: mcw64k3a2, hardcoded jumps from LCG state y, force a_1 = 0. \n";
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3a2();
   }
   tmp = clock() - tmp;
   printResultsTimex(n, jumpSize, tmp);
   // std::cout << "   MWC state = " << std::setw(4) << x0 << " "  << x1 << " " << x2 << " "  << x3 << "\n";

   std::cout << "\nGenerator name: mcw64k3a2s, hardcoded jumps from LCG state y + toMWC. \n";
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      LCGtoMWCHardk3a2s();
   }
   tmp = clock() - tmp;
   printResultsTimex(n, jumpSize, tmp);

   /*
   // This is a tentative implementation that uses GMP only, and not NTL.
   // mpz_t bmp, mmp, ymp, jumpMultmp, sigmamp, xx0, xx1, xx2, xx3;
   std::cout << "\n Generator name: mcw64k3a2, jumps from LCG state y + toMWC, using GMP. \n";
   std::cout << "\n  Time to make n = " << n << " jumps ahead, jump size = " << jumpSize << "\n";
   // getJumpAheadMult(jumpMultHard, b, m, jumpSize);
   jumpMultmp = conv<mpz_t>(jumpMultHard);
   xx0 = xx1 = xx2 = xx3 = 12345;
   for (int64_t j = 0; j <= 3; j++) xx[j] = 12345;
   MWCtoLCGStateLagk(yy, b, aa, xx);
   std::cout << "  state y = " << yy << "\n";
   ymp = conv<mpz_t>(yy);
   tmp = clock();
   for (int64_t i = 0; i < n; i++) {
      NTL::MulMod(yy, jumpMultHard, yy, m);
      ymp = conv<mpz_t>(yy);
      LCGtoMWCHardk3a2c();
   }
   tmp = clock() - tmp;
   std::cout << "  state y = " << ymp << "\n";
   // printResultsTime(xx, tmp);
   std::cout << "      " << std::setw(4) << xx0 << " "  << xx1 << " " << xx2 << " "  << xx3 << " "
        << std::fixed << std::setw(10) << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
  */

   return 0;

}

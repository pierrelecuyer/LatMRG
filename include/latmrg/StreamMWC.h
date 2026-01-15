/**
 * This does not belong to LatMRG.  ***
 * It is for testing multiple streams with a MWC generator.
 */

#ifndef STREAMMWC_H
#define STREAMMWC_H

#define MWC_U01_NOZERO

#include <string>
#include <iostream>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <stdint.h>
#include <math.h>
#include <iomanip> // For std::hex, std::setfill, std::setw
#include <algorithm> // Required for std::reverse

#define TYPES_CODE  ZD     // Int = ZZ, Real = double
#include <gmp.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"


namespace LatMRG {

using namespace std;

class StreamMWC {

public:

StreamMWC (const char *name = "");


StreamMWC (const uint64_t seed[4], const char *name = "");


static bool SetPackageSeed (const uint64_t seed[4]);


void ResetStartStream ();


void ResetStartSubstream ();


void ResetNextSubstream ();


bool SetSeed (const uint64_t seed[4]);


void AdvanceState (int64_t e, int64_t c);


void GetState (uint64_t seed[4]) const;


void WriteState () const;


void WriteStateFull () const;


uint64_t Rand64 ();


int64_t RandInt (int i, int j);


double RandU01 ();

// ------------------------------------

private:

uint64_t x1, x2, x3, c;
static __uint128_t tau, sum3;
// uint64_t sum = 0;

double Cg[4], Bg[4], Ig[4];

std::string name;

static uint64_t nextseed[4];
static const double twom53;
static const Int b;
static Int m;
static IntVec aa, xx;
static Int a0inv;
static Int yn;

// double U01 ();

};

//===========================================================================
// IMPLEMENTTION

// ***********************************************************************

// namespace {

const double twom53 = 1.0 / (double) ((uint64_t) 1 << 53);
const Int b = (NTL::power(Int(2), 64));


// From MWC1k3-0-60-58.res    a_0 = -1, a_1 = 0.
// uint64_t inline mwc64k3a2Xor() {
uint64_t inline StreamMWC::Rand64() {
   const uint64_t out = x1 ^ x3;
   tau = (0x2902ebc31ec3683 * (__uint128_t) x2 + 0x57baa090037 * (__uint128_t) x3) + c;
   x3 = x2;
   x2 = x1;
   x1 = tau;
   c = tau >> 64;
   return out;
}

//-------------------------------------------------------------------------
// Generate the next U01 random number.
//
double inline StreamMWC::RandU01 () {
   uint64_t block53 = Rand64() >> 11;
#ifdef MWC_U01_NOZERO
   if (block53 == 0) return RandU01();
#endif
   return block53 * twom53;
}


//*************************************************************************
// Public members of the class start here


//-------------------------------------------------------------------------
// The default seed of the package; will be the seed of the first
// declared StreamMWC, unless SetPackageSeed is called.
//
uint64_t StreamMWC::nextseed[4] = { 123456789, 123456789, 123456789, 123456789 };


// Transform LCG state `yy` to MWC state `xx`, for k=3, a_0=-1 and a_1=0.
static void inline LCGtoMWC() {
   Int sigma = yn;
   xx[0] = sigma  & (b-1);
   // x3 = conv<ulong>(sigma);
   sigma = (sigma - xx[0]) >> 64;
   xx[1] = (sigma) & (b-1);
   sigma = (sigma - xx[1]) >> 64;
   xx[2] = (sigma + a2 * xx[0]) & (b-1);
   xx[3] = (sigma - x2 + a2 * xx[0]) >> 64;
   nextseed[0] = conv<ulong>(xx[0]);
   nextseed[1] = conv<ulong>(xx[1]);
   nextseed[2] = conv<ulong>(xx[2]);
   nextseed[3] = conv<ulong>(xx[3]);
}

static void MWCtoLCG(Int &y, const IntVec &xx) {
   int64_t k = aa.length() - 1;
   Int bj = b;  // b^j
   Int dj;      // d_j
   y = -aa[0] * xx[0];  // = d_0
   for (int64_t j = 1; j < k; j++) {
      dj = Int(0);
      for (int64_t i = 0; i <= j; i++)
         dj -= aa[i] * xx[j-i];
      y += dj * bj;
      bj *= b;
   }
   y += xx[k] * bj;
}

//-------------------------------------------------------------------------
// constructor using jump ahead.
//
StreamMWC::StreamMWC (const char *s) : name (s) {
   for (int i = 0; i < 4; ++i) {
      Bg[i] = Cg[i] = Ig[i] = nextSeed[i];
   }
   LCGToMWC();
}


//-------------------------------------------------------------------------
// constructor with given seed.
//
StreamMWC::StreamMWC (const uint64_t seed[4], const char *s) : name (s) {
   for (int i = 0; i < 4; ++i) {
      Bg[i] = Cg[i] = Ig[i] = seed[i];
   }
}


//-------------------------------------------------------------------------
// Reset Stream to beginning of Stream.
//
void StreamMWC::ResetStartStream () {
   for (int i = 0; i < 4; ++i)
      Cg[i] = Bg[i] = Ig[i];
}


//-------------------------------------------------------------------------
// Reset Stream to beginning of SubStream.
//
void StreamMWC::ResetStartSubstream () {
   for (int i = 0; i < 4; ++i)
      Cg[i] = Bg[i];
}


//-------------------------------------------------------------------------
// Reset Stream to NextSubStream.
//
void StreamMWC::ResetNextSubstream () {
   for (int i = 0; i < 4; ++i)
       Cg[i] = Bg[i];
}


//-------------------------------------------------------------------------
void StreamMWC::SetPackageSeed (const uint64_t seed[4]) {
   for (int i = 0; i < 6; ++i)
      nextSeed[i] = seed[i];
}


//-------------------------------------------------------------------------
void StreamMWC::SetSeed (const uint64_t seed[4]) {
   for (int i = 0; i < 6; ++i)
      Cg[i] = Bg[i] = Ig[i] = seed[i];
}


//-------------------------------------------------------------------------
// if e > 0, let n = 2^e + c;
// if e < 0, let n = -2^(-e) + c;
// if e = 0, let n = c.
// Jump n steps forward if n > 0, backwards if n < 0.
//
void StreamMWC::AdvanceState (int64_t e, int64_t c)
{
    double B1[3][3], C1[3][3], B2[3][3], C2[3][3];

    if (e > 0) {
        MatTwoPowModM (A1p0, B1, m1, e);
        MatTwoPowModM (A2p0, B2, m2, e);
    } else if (e < 0) {
        MatTwoPowModM (InvA1, B1, m1, -e);
        MatTwoPowModM (InvA2, B2, m2, -e);
    }

    if (c >= 0) {
        MatPowModM (A1p0, C1, m1, c);
        MatPowModM (A2p0, C2, m2, c);
    } else {
        MatPowModM (InvA1, C1, m1, -c);
        MatPowModM (InvA2, C2, m2, -c);
    }

    if (e) {
        MatMatModM (B1, C1, C1, m1);
        MatMatModM (B2, C2, C2, m2);
    }

    MatVecModM (C1, Cg, Cg, m1);
    MatVecModM (C2, &Cg[3], &Cg[3], m2);
}


//-------------------------------------------------------------------------
void StreamMWC::GetState (uint64_t seed[4]) const {
   for (int i = 0; i < 4; ++i)
      seed[i] = static_cast<uint64_t> (Cg[i]);
}


//-------------------------------------------------------------------------
void StreamMWC::WriteState () const {
    cout << "The current state of the StreamMWC";
    if (name.size() > 0)
        cout << " " << name;
    cout << ":\n   Cg = { ";
    for (int i = 0; i < 3; i++) {
        cout << static_cast<uint64_t> (Cg [i]) << ", ";
    }
    cout << static_cast<uint64_t> (Cg [3]) << " }\n\n";
}

//-------------------------------------------------------------------------
// Generate the next random integer.
//
int64_t StreamMWC::RandInt (int64_t low, int64_t high) {
    return low + static_cast<int64_t> ((high - low + 1.0) * RandU01 ());
};

// } // end of anonymous namespace

} // End namespace

#endif






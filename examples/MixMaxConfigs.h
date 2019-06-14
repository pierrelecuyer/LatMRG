#ifndef LATMRG_MIXMAXCONFIGS_H
#define LATMRG_MIXMAXCONFIGS_H

#include "NTL/ZZ.h"

#include "latticetester/ntlwrap.h"

/**
 * This is a header file containing constants and definitions for the MixMax
 * example.
 * */

typedef NTL::ZZ Int;
typedef NTL::vector<NTL::ZZ> IntVec;
typedef double Dbl;

struct MixMaxConfig {
  Int m, d, c;

  int k;

  MixMaxConfig(Int m, int k, Int d, Int c) {
    this->m = m;
    this->k = k;
    this->d = d;
    this->c = c;
  }
};

const int num_confs = 3;
MixMaxConfig Configs[3] = {MixMaxConfig(NTL::power_ZZ(2,61)-1, 8, Int(0), NTL::power_ZZ(2,53)+1),
MixMaxConfig(NTL::power_ZZ(2,61)-1, 17, Int(0), NTL::power_ZZ(2,36)+1),
MixMaxConfig(NTL::power_ZZ(2,61)-1, 240, Int(487013230256099140), NTL::power_ZZ(2,51)+1)};

#endif

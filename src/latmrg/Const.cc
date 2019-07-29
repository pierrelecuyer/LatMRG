// This file is part of LatMRG.
//
// LatMRG
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string>
#include <iostream>

#include "latmrg/Const.h"


using namespace std;

namespace LatMRG
{

  //===========================================================================

  string toStringGen (GenType gen)
  {
    switch (gen) {
      case MWC:
        return "MWC";
      case LCG:
        return "LCG";
      case MRG:
        return "MRG";
      case MMRG:
        return "MMRG";
      case KOROBOV:
        return "KOROBOV";
      case RANK1:
        return "RANK1";
      case COMBO:
        return "COMBO";
      default:
        return "***** GenType: IMPOSSIBLE CASE ";
    }
  }

  //============================================================================

  int toGenString(GenType& type, const string& type_str) {
    if (type_str == "MRG") {
      type = MRG;
    } else if (type_str == "LCG") {
      type = LCG;
    } else if (type_str == "MWC") {
      type = MWC;
    } else if (type_str == "KOROBOV") {
      type = KOROBOV;
    } else if (type_str == "RANK1") {
      type = RANK1;
    } else if (type_str == "MMRG") {
      type = MMRG;
    } else if (type_str == "COMBO") {
      type = COMBO;
    } else {
      std::cerr << type_str << " is not a GenType.\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  string toStringLattice (LatticeType lat)
  {
    switch (lat) {
      case FULL:
        return "FULL";
      case RECURRENT:
        return "RECURRENT";
      case ORBIT:
        return "ORBIT";
      case PRIMEPOWER:
        return "PRIMEPOWER";
      default:
        return "***** LatticeType: IMPOSSIBLE CASE ";
    }
  }

  //===========================================================================

  string toStringDecomp (DecompType deco)
  {
    switch (deco) {
      case DECOMP:
        return "DECOMP";
      case DECOMP_WRITE:
        return "DECOMP_WRITE";
      case DECOMP_READ:
        return "DECOMP_READ";
      case DECOMP_PRIME:
        return "DECOMP_PRIME";
      default:
        return "***** DecompType: IMPOSSIBLE CASE ";
    }
  }

  //===========================================================================

  string toStringImplemCond (ImplemCond cond)
  {
    switch (cond) {
      case NO_COND:
        return "NO_COND";
      case APP_FACT:
        return "APP_FACT";
      case EQUAL_COEF:
        return "EQUAL_COEF";
      case ZERO_COEF:
        return "ZERO_COEF";
      case POWER_TWO:
        return "POWER_TWO";
      default:
        return "***** ImplemCond: IMPOSSIBLE CASE ";
    }
  }


  //===========================================================================

  string toStringSearchMethod (SearchMethod method)
  {
    switch (method) {
      case EXHAUST:
        return "EXHAUST";
      case RANDOM:
        return "RANDOM";
      default:
        return "***** SearchMethod: IMPOSSIBLE CASE ";
    }
  }

  //===========================================================================

  int toCriterionString(LatticeTester::CriterionType& criter, const std::string& crit_str) {
    if (crit_str == "LENGTH") {
      criter = LatticeTester::LENGTH;
    } else if (crit_str == "SPECTRAL") {
      criter = LatticeTester::SPECTRAL;
    } else if (crit_str == "BEYER") {
      criter = LatticeTester::BEYER;
    } else if (crit_str == "PALPHA") {
      criter = LatticeTester::PALPHA;
    } else if (crit_str == "BOUND_JS") {
      criter = LatticeTester::BOUND_JS;
    } else {
      std::cerr << crit_str << " is not a CriterionType.\n";
      return 1;
    }
    return 0;
  }

  //============================================================================

  int toRedString(LatticeTester::PreReductionType& red, const std::string& red_str) {
    if (red_str == "FULL") red = LatticeTester::FULL;
    else if (red_str == "BKZ") red = LatticeTester::BKZ;
    else if (red_str == "DIETER") red = LatticeTester::DIETER;
    else if (red_str == "LLL") red = LatticeTester::LLL;
    else if (red_str == "NOPRERED") red = LatticeTester::NOPRERED;
    else {
      std::cerr << red_str << " is not a PreReductionType\n";
      return 1;
    }
    return 0;
  }

  //===========================================================================

  int toNormaString(LatticeTester::NormaType& norma, const std::string& norma_str) {
    if (norma_str == "BESTLAT") {
      norma = LatticeTester::BESTLAT;
    } else if (norma_str == "BESTBOUND") {
      norma = LatticeTester::BESTBOUND;
    } else if (norma_str == "LAMINATED") {
      norma = LatticeTester::LAMINATED;
    } else if (norma_str == "ROGERS") {
      norma = LatticeTester::ROGERS;
    } else if (norma_str == "MINK") {
      norma = LatticeTester::MINK;
    } else if (norma_str == "MINKL1") {
      norma = LatticeTester::MINKL1;
    } else if (norma_str == "L1") {
      norma = LatticeTester::L1;
    } else if (norma_str == "L2") {
      norma = LatticeTester::L2;
    } else if (norma_str == "NONE") {
      norma = LatticeTester::NONE;
    } else {
      std::cerr << norma_str << " is not a NormaType\n";
      return 1;
    }

    return 0;
  }

}

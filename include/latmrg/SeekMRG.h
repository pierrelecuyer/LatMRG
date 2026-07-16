// This file is part of LatMRG.
//
// Copyright (C) 2012-2022  The LatMRG authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit? de Montr?al.
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

#ifndef LATMRG_SEEKMRG_H
#define LATMRG_SEEKMRG_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>

#include "latmrg/PrimesFinder.h"
#include "latmrg/ConfigSeek.h"
#include "latmrg/FigureOfMeritData.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsUniform.h"

#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"
#include "latmrg/MWCLattice.h"

namespace LatMRG {

    template<typename Int, typename Real>
    class SeekMRG: public Seek<MRGLattice<Int, Real>> {

        public:

        using Seek<MRGLattice<Int, Real>>::Seek;

        MRGLattice<Int, Real>* nextGenerator() override;

        MRGLattice<Int, Real>* nextGeneratorRandom() override;

        private:

    }
  }

    //============================================================================
// Implementation


  /**
  * This method yields the next generator of an MRG lattice based on the chosen 
  * configuration.
  */
  MRGLattice<Int, Real>* SeekMRG::nextGenerator()
  {
    auto* comp = conf.genComponents[0];

    Int range, val;
    Int tmp = NTL::to_ZZ(currentGen);

    int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    // decode currentGen into multi-index
    for (int i = 1; i <= k; i++) {
      range = comp->getHighBoundary(i) - comp->getLowBoundary(i) + 1;
      val = tmp % range;
      tmp /= range;
      aa[i] = comp->getLowBoundary(i) + val;
    }
    
    if (currentGen < conf.max_gen && currentGen < comp->getNoMultipliers())
    {
      ++currentGen;
      return new MRGLattice<Int, Real>(comp->getModulus(), aa, conf.maxdim);
    }
    // Otherwise return null pointer
    return nullptr;
  }

   /**
  * This method yields a random generator of an MRG lattice within the range
  * defined by the chosen configuration.
  */  
  MRGLattice<Int, Real>* SeekMRG::nextGeneratorRandom()
  {
    const auto* comp = conf.genComponents[0];
    const int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    for (int i = 1; i <= k; i++) {
      // CW: There seems to be some problem with RandInt after a few calls. Need to check this at a later point.
      aa[i] = randInt(comp->getLowBoundary(i), comp->getHighBoundary(i));
    }

    if (currentGen < conf.max_gen)
    {
      ++currentGen;
      return new MRGLattice<Int, Real>(comp->getModulus(), aa, conf.maxdim);  
    }
  }

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

#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

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

  /**
 * \class Seek
 *
 * This class performs an exhaustive or random search for linear
 * recurrence generators satisfying prescribed constraints and ranks
 * the resulting generators according to a chosen figure of merit.
 *
 * A search consists of generating candidate lattices one by one,
 * optionally verifying that the corresponding generator has maximal
 * period, and computing a figure of merit for each admissible lattice.
 * Depending on the configuration, the figure of merit is evaluated
 * either for the primal or the dual lattice. The best generators
 * encountered during the search are retained throughout the execution.
 *
 * The constructor initializes the search from a ConfigSeek object.
 * In particular, it stores the search configuration, constructs the
 * required figure of merit objects, and initializes the auxiliary
 * objects needed during the search.
 *
 * The method performSeek() is the main facility of the
 * class. It repeatedly generates candidate lattices, evaluates their
 * figure of merit, and keeps the best generators found during the
 * search. Note that the candidate lattice depend on the configuration
 * as well as on the lattice type under consideration. The search terminates 
 * when no further admissible generators are available or when the 
 * prescribed CPU time limit is reached
 *
 * Candidate generators are produced by the member functions
 * \nextGenerator() and nextGeneratorRandom(). These functions
 * generate the next admissible generator according to the search
 * configuration, either by exhaustive enumeration of the prescribed
 * search space or by random sampling from the admissible parameter
 * ranges. The implementation of these functions is specialized for the
 * supported lattice types [CW: In the subclasses].
 */

  template<typename Lat> class Seek {
        
    private:
        
    /**
    * The current lattice which is analyzed.
    */
    std::unique_ptr<Lat> lat;

    /**
    * The configuration for the seek.
    */
    ConfigSeek<Int, Real> conf;

    /**
    * Figure of merit objects for the primal and dual lattice.
    */        
    FigureOfMeritM<Int, Real> fomPrimal;
    FigureOfMeritDualM<Int, Real> fomDual;

    /**
     * Object for creating an MRG lattice
    */
    MRGComponent<Int> mrg;

    /**
    * Program time
    */
    Chrono timer; 

    /**
    * Counter for the current RNG.
    */
    int currentGen = 0;       
         
    public:

    /**
    * Initialization of objects is based on the ConfigSeek parameter.
    */
    Seek(const ConfigSeek<Int, Real>& config) : conf(config) 
                                          , fomPrimal(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst)
                                          , fomDual(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst) 
                                          , mrg(config.genComponents[0]->getModulus(), config.genComponents[0]->getOrder())
                                          {                                             
                                          }
    
                                          
    /**
    * This method looks for the next generator based on the chosen configuration.
    * It uses the variable currentGen to calculate the next generator in the list.
    * nextGenerator() is performing an exhaustive search. This method needs to
    * be implemented separately for the different types of lattices.
    */
    Lat* nextGenerator();

    /**
    * As nextGenerator() but it chooses a random next generator which is within
    * the range defined by the configuration.
    */
    Lat* nextGeneratorRandom();


    /**
    * This method performs the seek based on the current configuration.
    * The parameter *generator defines the method of the seek, e.g.,
    * a random or an exhaustive search.
    */
    int performSeek(Lat* (Seek::*generator)());
    
    /**
    * This method is supposed to report the progress of the current  
    * search. It is still the version of the old LatMRG and currently
    * not working.
    */
    int print_progress(int old);

  };

    
//============================================================================
// Implementation


  /**
  * This method yields the next generator of an MRG lattice based on the chosen 
  * configuration.
  */
  template<> MRGLattice<Int, Real>* Seek<MRGLattice<Int, Real>>::nextGenerator()
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
   * Same for MWC Lattices
   */
  template<> MWCLattice<Int, Real>* Seek<MWCLattice<Int, Real>>::nextGenerator()
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

    Int b = NTL::power(Int(2), comp->getPowMod());

    if (currentGen < conf.max_gen && currentGen < comp->getNoMultipliers())
    {
      ++currentGen;
      return new MWCLattice<Int, Real>(b, aa, conf.maxdim, conf.maxdim, conf.maxdim);
    }

    return nullptr;
  }

  /**
  * This method yields a random generator of an MRG lattice within the range
  * defined by the chosen configuration.
  */  
  template<> MRGLattice<Int, Real>* Seek<MRGLattice<Int, Real>>::nextGeneratorRandom()
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

    return nullptr;
  }
  
  /**
   * This method yields a random MWC Lattice currently, it is restricting the number of bits used
   * per component of the multiplier. However, there are no restrictions on the range of the 
   * coefficients.
   */
  template<> MWCLattice<Int, Real>* Seek<MWCLattice<Int, Real>>::nextGeneratorRandom()
  {
    auto* comp = conf.genComponents[0];
    int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);
    // Initialized
    aa[0] = -1;
    for (int64_t j = 1; j < k; j++)
      aa[j] = 0;
    Int b = NTL::power(Int(2), comp->getPowMod());

    // Get random MWC generator
    if (comp->getRandomBits(0) > 0) {
         aa[0] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(0)));
         if ((aa[0] % 2) == 0) aa[0] += 1;  // a_0 must be odd.
      }
    aa[k] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(k)));  // `ek` random bits for a_k.
    for (int64_t j = 1; j < min(asMWC(comp)->numaj, k); j++)
        aa[k-j] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(j)));  // `ej` random bits for a_{k-j}.
    for (int64_t j = 1; j < -asMWC(comp)->numaj; j++)  // To impose equal coefficients, for testing.
        aa[k-j] = aa[k];
    
    if (currentGen < conf.max_gen)
    {
      ++currentGen;
      return new MWCLattice<Int, Real>(b, aa, conf.maxdim, conf.maxdim, conf.maxdim);
    }

    return nullptr;
  }
  
  /**
  * This method is supposed to report the progress of the current  
  * search. It is still the version of the old LatMRG and currently
  * not working.
  */
  template<typename Lat> int Seek<Lat>:: print_progress(int old) {    
    int per_80 = 80 * timer.val(Chrono::SEC) / conf.timeLimit;
    if (per_80 > 80) per_80 = 80;
    if (per_80 < 0) per_80 = 0;
    // We do not print for no reason as this slows the program a lot.
    if (per_80 <= old) return old;
    std::cout << "Program progress: [";
    for (int i = 0; i < per_80; i++) std::cout << "#";
    for (int i = per_80; i < 80; i++) std::cout << " ";
    std::cout << "] ";
    std::cout << std::setw(3) << int(per_80/80.0*100) << " %\r" << std::flush;
    return per_80;
  }

  /**
  * This method performs the seek based on the current configuration.
  * The parameter *generator defines the method of the seek, e.g.,
  * a random or an exhaustive search.
  */
  template<typename Lat> int Seek<Lat>::performSeek(Lat* (Seek::*generator)())  {  
    assert(!conf.genComponents.empty());  
    const auto* comp = conf.genComponents[0];
    bool maxper = true;
    int old = 0;
    // Launching the tests
    if (conf.progress) {
      old = print_progress(-1);
    }
    MeritList<Lat> bestLattice(conf.configFOM.no_bestGen, conf.configFOM.best);
    timer.init();
    
    IntLattice<Int, Real> proj(comp->getModulus(), conf.configFOM.t.length(), conf.configFOM.norm);
    
    FigureOfMeritData<Lat> fomData;

    // Preparation for being able to check primitivity
    setModulusIntP<Int>(comp->getModulus());
    
    // Loop through all lattices
    do {      
      lat.reset((this->*generator)());
      if (!lat) continue;   
      // If the period must be maximal, test if the specific component has max period.
      // Otherwise continue.
      // This needs to be amended if we want to organize the code by subclasses for the specific RNGs.
      if (comp->onlyMaxPeriod())
      {
        // Check the different cri
        if (conf.genType == MRG) {
          maxper = mrg.maxPeriod(lat->getaa());
        }
        else if (conf.genType == MWC) {
            maxper = mIsSafePrime(lat->getModulus(), 100) && maxPeriodHalfMWC (lat->getModulus(), comp->getPowMod());
        }
        if (!maxper)
            continue;
      }
      
      // Variable which depends on whether the FoM is calculated for the primal / dual lattice.
      // Avoids duplicate code.
      auto& fom = conf.configFOM.dualLattice ? fomDual : fomPrimal;     
      fomData.setMerit(fom.computeMerit(*lat, proj));
      fomData.setLattice(lat.get());
      fomData.setMeritProj(fom.getMinMeritProj());
      fomData.setMeritSqlen(fom.getMinMeritSqlen());

      bestLattice.add(fomData); 
      // Update the lower bound for the FoM to equal the smallest stored FoM
      if (conf.configFOM.best)
        fom.setLowBound(bestLattice.getSmallestMerit());

      conf.configFOM.num_gen++;
      conf.configFOM.currentMerit = bestLattice.getMerit();
      if (conf.progress) old = print_progress(old);
    } while (!timer.timeOver(conf.timeLimit) && lat);
     
    return 0;
  }

}// End namespace LatMRG
#endif

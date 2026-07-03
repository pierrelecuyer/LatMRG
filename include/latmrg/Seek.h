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

namespace LatMRG {

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
                                          {                                             
                                          }
    
                                          
    /**
    * This method looks for the next generator based on the chosen configuration.
    * It uses the variable currentGen to calculate the next generator in the list.
    * nextGenerator() is performing an exhaustive search. This method needs to
    * be implemented separately for the different types of lattices.
    */
    MRGLattice<Int, Real>* nextGenerator();

    /**
    * As nextGenerator() but it chooses a random next generator which is within
    * the range defined by the configuration.
    */
    MRGLattice<Int, Real>* nextGeneratorRandom();

    /**
    * This method performs the seek based on the current configuration.
    * The parameter *generator defines the method of the seek, e.g.,
    * a random or an exhaustive search.
    */
    int performSeek(MRGLattice<Int, Real>* (Seek::*generator)());  
    
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
  template<typename Lat> MRGLattice<Int, Real>* Seek<Lat>::nextGenerator()
  {
    auto* comp = conf.genComponents[0];

    Int range, val;
    Int tmp = NTL::to_ZZ(currentGen);

    int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> a;
    a.SetLength(k + 1);

    // decode currentGen into multi-index
    for (int i = 1; i <= k; i++) {
      range = comp->getHighBoundary(i) - comp->getLowBoundary(i) + 1;
      val = tmp % range;
      tmp /= range;
      a[i] = comp->getLowBoundary(i) + val;
    }
    
    if (currentGen < conf.configFOM.max_gen && currentGen < comp->getNoMultipliers())
    {
      ++currentGen;
      return new MRGLattice<Int, Real>(comp->getModulus(), a, conf.maxdim);
    }
    // Otherwise return null pointer
    return nullptr;
  }

  /**
  * This method yields a random generator of an MRG lattice within the range
  * defined by the chosen configuration.
  */  
  template<typename Lat> MRGLattice<Int, Real>* Seek<Lat>::nextGeneratorRandom()
  {
    const auto* comp = conf.genComponents[0];
    const int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> a;
    a.SetLength(k + 1);

    for (int i = 1; i <= k; i++) {
      // CW: There seems to be some problem with RandInt after a few calls. Need to check this at a later point.
      a[i] = randInt(comp->getLowBoundary(i), comp->getHighBoundary(i));
    }

    if (currentGen < conf.configFOM.max_gen)
    {
      currentGen++;
      return new MRGLattice<Int, Real>(comp->getModulus(), a, conf.maxdim);
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
  template<typename Lat> int Seek<Lat>::performSeek(MRGLattice<Int, Real>* (Seek::*generator)())  {  
    assert(!conf.genComponents.empty());  
    const auto* comp = asMRG(conf.genComponents[0]);
    int old = 0;
    // Launching the tests
    if (conf.progress) {
      old = print_progress(-1);
    }
    MeritList<Lat> bestLattice(conf.configFOM.max_gen, conf.configFOM.best);
    timer.init();
    
    IntLattice<Int, Real> proj(comp->getModulus(), conf.configFOM.t.length(), conf.configFOM.norm);
    
    FigureOfMeritData<Lat> fomData;

    // Preparation for being able to check primitivity
    setModulusIntP<Int>(comp->getModulus());
    MRGComponent<Int> mrg(comp->getModulus(), comp->order);
    
    // Loop through all MRGs
    do {      
      lat.reset((this->*generator)());
      if (!lat) continue;   
      // If the period must be maximal, test if the specific component has max period.
      // Otherwise continue.
      if (comp->permaxPrime)
      {
        if (!mrg.maxPeriod(lat->getaa()))
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
      conf.configFOM.num_gen++;
      conf.configFOM.currentMerit = bestLattice.getMerit();
      if (conf.progress) old = print_progress(old);
    } while (!timer.timeOver(conf.timeLimit) && lat);
     
    return 0;
  }

}// End namespace LatMRG
#endif

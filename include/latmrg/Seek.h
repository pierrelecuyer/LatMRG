#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

#include <cassert>

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
    Lat* lat = 0;

    /**
    * The configuration of the seek
    */
    ConfigSeek<Int, Real> conf;

    /**
    * Firgure of merit object for the primal and dual lattice.
    */        
    FigureOfMeritM<Int, Real> fomPrimal;
    FigureOfMeritDualM<Int, Real> fomDual;

    /**
    * Program time
    */
    Chrono timer; 

    /**
    * Counter for the current RNG
    */
    Int currentGen;       
         
    public:

    /**
    * Initializition of objects is based on the ConfigSeek parameter 
    */
    Seek(ConfigSeek<Int, Real>& config) : conf(config), 
                                          fomPrimal(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst),
                                          fomDual(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst) { };
    
                                          
    /**
    * This method looks for the next generator based on the chosen configuration.
    * It uses the variable currentGen to calculate the next generator in the list.
    * nextGenerator() is performing an exhaustive search. This method needs to
    * be implemented separately for the different types of lattices.
    */
    template<typename Int, typename Real> MRGLattice<Int, Real>* nextGenerator();

    /**
    * As nextGenerator() but it chooses a random next generator which is within
    * the range defned by the configuration
    */
    template<typename Int, typename Real> MRGLattice<Int, Real>* nextGeneratorRandom();

    /**
    * This method performs the seek based on the current configuration.
    * The parameter *generator defines the method of the seek, e.g.,
    * a random or an exhaustive search
    */
    int PerformSeek(MRGLattice<Int, Real>* (Seek::*generator)());  
    
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
  template<typename Lat> template<typename Int, typename Real> MRGLattice<Int, Real>*Seek<Lat>::nextGenerator()
  {
    auto* comp = asMRG(this->conf.genComponents[0]);

    Int range, val;
    Int tmp = currentGen;

    int k = comp->order;
    NTL::Vec<NTL::ZZ> a;
    a.SetLength(k + 1);

    // decode currentGen into multi-index
    for (int i = 1; i <= k; i++) {
      range = comp->highBoundaries[i] - comp->lowBoundaries[i] + 1;
      Int val = tmp % range;
      tmp /= range;
      a[i] = comp->lowBoundaries[i] + val;
    }
    
    if (currentGen < this->conf.configFOM.max_gen && currentGen < comp->getNoMultipliers())
    {
      currentGen++;
      return new MRGLattice<Int, Real>(comp->modulus, a, this->conf.maxdim);
    }
    // Otherwise return null pointer
    return nullptr;
  }

  /**
  * This method yields a random generator of an MRG lattice within the range
  * defined by the chosen configuration.
  */  
  template<typename Lat> template<typename Int, typename Real> MRGLattice<Int, Real>*Seek<Lat>::nextGeneratorRandom()
  {
    auto* comp = asMRG(this->conf.genComponents[0]);
    int k = comp->order;
    NTL::Vec<NTL::ZZ> a;
    a.SetLength(k + 1);

    for (int i = 1; i <= k; i++) {
      // CW: There seems to be some problem with RandInt after a few calls. Need to check this at a later point.
      a[i] = randInt(comp->lowBoundaries[i], comp->highBoundaries[i]);
    }

    if (currentGen < this->conf.configFOM.max_gen)
    {
      currentGen++;
      std::cout << currentGen << "\n";
      return new MRGLattice<Int, Real>(comp->modulus, a, this->conf.maxdim);
    }

    return nullptr;
  }
  
  /**
  * This method is supposed to report the progress of the current  
  * search. It is still the version of the old LatMRG and currently
  * not working.
  */
  template<typename Lat> int Seek<Lat>:: print_progress(int old) {    
    int per_80 = timer.val(Chrono::SEC)/this->conf.timeLimit * 80;
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
  * a random or an exhaustive search
  */
  template<typename Lat> int Seek<Lat>::PerformSeek(MRGLattice<Int, Real>* (Seek::*generator)())  {    
    int old = 0;
    // Launching the tests
    if (this->conf.progress) {
      old = print_progress(-1);
    }
    MeritList<Lat> bestLattice(this->conf.configFOM.max_gen, this->conf.configFOM.best);
    timer.init();
    
    IntLattice<Int, Real> proj(this->conf.genComponents[0]->getModulus(), this->conf.configFOM.t.length(), this->conf.configFOM.norm);
    IntLattice<Int, Real> m_lattice(this->conf.genComponents[0]->getModulus(), this->conf.configFOM.t.length(), this->conf.configFOM.norm);
    
    FigureOfMeritData<Lat> fomData;

    // Preparation for being able to check primitivity
    setModulusIntP<Int>(this->conf.genComponents[0]->getModulus());
    MRGComponent<Int> mrg(this->conf.genComponents[0]->getModulus(), asMRG(this->conf.genComponents[0])->order);
    
    // Loop through all MRGs
    do {      
      if (lat != NULL) delete lat;
      lat = (this->*generator)();
      if (lat == NULL) continue;   
      // If the period must be maximal, test if the specific component has max period.
      // Otherwise continue.
      if (asMRG(this->conf.genComponents[0])->permaxPrime)
      {
        if (!mrg.maxPeriod(lat->getaa()))
            continue;
      }
      // Calculate the FOM for the primal / dual
      if (this->conf.configFOM.dualLattice) {        
        fomData.setMerit(this->fomDual.computeMerit(*lat, proj));
        fomData.setLattice(lat);
        fomData.setMeritProj(this->fomDual.getMinMeritProj());
        fomData.setMeritSqlen(this->fomDual.getMinMeritSqlen());;
      } 
      else {
        fomData.setMerit(this->fomPrimal.computeMerit(*lat, proj));
        fomData.setLattice(lat);
        fomData.setMeritProj(this->fomPrimal.getMinMeritProj());
        fomData.setMeritSqlen(this->fomPrimal.getMinMeritSqlen());
      }
      bestLattice.add(fomData); 
      this->conf.configFOM.num_gen++;
      this->conf.configFOM.currentMerit = bestLattice.getMerit(); 
      if (this->conf.progress) old = print_progress(old);
    } while (!timer.timeOver(this->conf.timeLimit) && this->lat);
     
    return 0;
  }

}// End namespace LatMRG
#endif

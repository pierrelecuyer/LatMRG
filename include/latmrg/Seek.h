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

        Lat* lat = 0;
        ConfigSeek<Int, Real> conf;
        ReducerBB<Int, Real> red;
        WeightsUniform weights;
        FigureOfMeritM<Int, Real> fomPrimal;
        FigureOfMeritDualM<Int, Real> fomDual;
        Chrono timer; // program timer

        Int currentGen; // CW: counter to get the nummber of the current generator
        
         
        public:

        Seek(ConfigSeek<Int, Real>& conf) :     conf(conf), 
                                                red(conf.maxdim), 
                                                weights(1.0),
                                                fomPrimal(conf.configFOM.t, weights, *conf.configFOM.norma, &red, true),
                                                fomDual(conf.configFOM.t, weights, *conf.configFOM.norma, &red, true) { };

        int PerformSeek();  
        
        template<typename Int, typename Real> MRGLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf);

        int print_progress(int old);

    };

    
//============================================================================
// Implementation

  template<typename Lat> template<typename Int, typename Real> MRGLattice<Int, Real>*Seek<Lat>::nextGenerator(ConfigSeek<Int, Real>& conf)
  {
     auto* comp = asMRG(conf.genComponents[0]);

      int k = comp->order;
      NTL::Vec<NTL::ZZ> a;
      a.SetLength(k + 1);

      Int tmp = currentGen;

      // decode currentGen into multi-index
      for (int i = 1; i <= k; i++) {
        Int lo = comp->lowBoundaries[i];
        Int hi = comp->highBoundaries[i];
        Int range = hi - lo + 1;
        Int val = tmp % range;
        tmp /= range;
        a[i] = lo + val;
      }
      std::cout << "a: " << a << "\n";
      // compute total number of combinations
      Int total;
      total = 1;
      for (int i = 1; i <= k; i++) {
          total *= (comp->highBoundaries[i] - comp->lowBoundaries[i] + 1);
      }

      if (currentGen < conf.configFOM.max_gen && currentGen < total)
      {
        currentGen++;
        return new MRGLattice<Int, Real>(comp->modulus, a, conf.maxdim);
      }

      return nullptr;
    /*
    NTL::Vec<NTL::ZZ> a;
    a.SetLength(asMRG(conf.genComponents[0])->order+1);
    if (currentGen < this->conf.configFOM.max_gen && asMRG(conf.genComponents[0])->lowBoundaries[1] + currentGen <= asMRG(conf.genComponents[0])->highBoundaries[1] )
      {
        a[1] = asMRG(conf.genComponents[0])->lowBoundaries[1] + currentGen;
        currentGen++;
        return new MRGLattice<Int, Real>(asMRG(conf.genComponents[0])->modulus, a, conf.maxdim);
      }
    else return nullptr;
    */
  }

  template<typename Lat> int Seek<Lat>:: print_progress(int old) {
    
    int per_80 = timer.val(Chrono::SEC)/conf.timeLimit * 80;
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

template<typename Lat> int Seek<Lat>::PerformSeek()  {
    
    int old = 0;
    // Launching the tests
    if (conf.progress) {
      old = print_progress(-1);
    }
    MeritList<Lat> bestLattice(this->conf.configFOM.max_gen, this->conf.configFOM.best);
    timer.init();

    Int modulus; // CW: interim: to be replaced later on by a variable stored in conf
    
    IntLattice<Int, Real> proj(modulus, this->conf.configFOM.t.length(), this->conf.configFOM.norm);
    
    MRGLattice<Int, Real> m_lattice(modulus, 1, this->conf.maxdim);

    FigureOfMeritData<Lat> fomData;
    
    do {      
      if (lat != NULL) delete lat;
      lat = nextGenerator(conf);
      if (lat == NULL) continue;           
      if (conf.configFOM.dualLattice) {        
        fomData.setMerit(this->fomDual.computeMerit(*lat, proj));
        std::cout << old << ": " << fomData.getMerit() << "\n";
        fomData.setLattice(lat);
        fomData.setMeritProj(this->fomDual.getMinMeritProj());
        fomData.setMeritSqlen(this->fomDual.getMinMeritSqlen());;
      } 
      else {
        fomData.setMerit(this->fomPrimal.computeMerit(*lat, proj));
        std::cout << old << ": " << fomData.getMerit() << "\n";
        fomData.setLattice(lat);
        fomData.setMeritProj(this->fomPrimal.getMinMeritProj());
        fomData.setMeritSqlen(this->fomPrimal.getMinMeritSqlen());
      }
      bestLattice.add(fomData);
      conf.configFOM.num_gen++;
      conf.configFOM.currentMerit = bestLattice.getMerit(); 
      if (conf.progress) old = print_progress(old);
   } while (!timer.timeOver(conf.timeLimit) && this->lat);
     
    return 0;
}

}// End namespace LatMRG
#endif

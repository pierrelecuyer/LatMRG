#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

#include <cassert>

#include "latmrg/PrimesFinder.h"
#include "latmrg/ConfigSeek.h"
#include "latmrg/Projections.h"
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
    MeritList<Lat> bestLattice(this.conf.configFOM.max_gen, this.conf.configFOM.best);
    timer.init();

    Int modulus; // CW: interim: to be replaced later on by a variable stored in conf

    FigureOfMeritData<Lat> fomData;
    IntLattice<Int, Real> proj(modulus, conf.configFOM.t.length(), conf.configFOM.norm);
    
    do {
      if (lat != NULL) delete lat;
      lat = nextGenerator(conf);
      if (lat == NULL) continue;
      if (conf.configFOM.dualLattice) {
        fomData.getMerit() = this->fomDual.computeMerit(lat, proj);
        fomData.getLattice() = lat;
        fomData.getMeritProj() = this->fomDual.getMinMeritProj();
        fomData.getMeritSqlen() = this->fomDual.getMinMeritSqlen();
      } 
      else {
        fomData.getMerit() = this->fomPrimal.computeMerit(lat, proj);
        fomData.getLattice() = lat;
        fomData.getMeritProj() = this->fomPrimal.getMinMeritProj();
        fomData.getMeritSqlen() = this->fomPrimal.getMinMeritSqlen();
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

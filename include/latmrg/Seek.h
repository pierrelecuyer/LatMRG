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

        
        ConfigSeek<Int, Real> conf;
        Int modulus = Int(1021); 
        ReducerBB<Int, Real> red;
        WeightsUniform weights;
        FigureOfMeritM<Int, Real> fomPrimal;
        FigureOfMeritDualM<Int, Real> fomDual;

         
        public:

        Seek(ConfigSeek<Int, Real>& conf) :     conf(conf), 
                                                red(conf.maxdim), 
                                                weights(1.0),
                                                normLattice(log(modulus), 1, conf.maxdim, conf.configFOM.norm),
                                                fomPrimal(conf.configFOM.t, weights, *conf.configFOM.norma, &red, true),
                                                fomDual(conf.configFOM.t, weights, *conf.configFOM.norma, &red, true) { };// ,
                                                // fomPrimal(t, weights, norm, &red, true) { };

        int PerformSeek();  

    };

    
//============================================================================
// Implementation

template<typename Lat> int Seek<Lat>::PerformSeek()  {
    MeritList<Lat> bestLattice(this.conf.configFOM.max_gen, this.conf.configFOM.best);
    return 0;
}

}// End namespace LatMRG
#endif

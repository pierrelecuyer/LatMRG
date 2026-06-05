#ifndef LATMRG_CONFIGSEEK_H
#define LATMRG_CONFIGSEEK_H

#include <ctime>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"

#include "latmrg/EnumTypes.h"
#include "latmrg/MWCLattice.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"
#include "latmrg/Projections.h"

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace std;
using namespace LatMRG;
using namespace LatticeTester::Random;


/**
 * Abstract structure to store the information required for a component in a seek.
 */
template<typename Int, typename Real> struct ConfigSeekComponent
{
    virtual ~ConfigSeekComponent() = default;

    std::int64_t nodesBB = 100000000000;

    bool use_dual = true;

    // Projection
    Projections* proj = nullptr;
    // Common component data
    Int b;

    // MRGComponent stores the stuff we might want to know, such as the modulo
    // the order and even the coefficients
    int num_comp = 1;
    std::vector<std::string> search_mode;
    //std::vector<bool> period; // CW: was a vector - I do not think this is true since we the ConfigSeekComponent contains the individual components
    bool period; // CW: new
    std::vector<MRGComponent<Int>*> fact;

    /**
     * Factory function.
     * Creates the correct subclass automatically.
     */
    static ConfigSeekComponent<Int, Real>* create(GenType type);
};


/**
 * Configuration for a MRG component.
 */
template<typename Int, typename Real> struct ConfigSeekMRG : ConfigSeekComponent<Int, Real>
{
    Int modulus;   // The modulus m = b^e + r
    Int base;      // The base b
    Int exponent;  // The exponent e
    Int rest;      // The rest r
    long order;    // The order k
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j 
    IntVec highBoundaries; // Default values are 1 and m-1 
    bool appFact;   // true iff we enforce the approximate factoring condition (default = false)
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
    bool permaxPrime;  // When true, m must be prime and we retain only max-period MRGs (default = false)
    DecompType howFactorh; // Indicates how to factorize h=(m-1)/2
    IntVec factorsh;      // The factors of h, repeated according to multiplicity, when known
    DecompType howFactor;  // Indicates how to factorize r=(m^k-1)/(m-1)
    IntVec factorsr;      // The factors of r, repeated according to multiplicity, when known
    bool permaxPow2;   // True iff m is a power of 2, b=2, r=0, k=1, and we want max period for the LCG
                       // Default is false
};


/**
 * Configuration for a MWC component.
 */
template<typename Int, typename Real>
struct ConfigSeekMWC : ConfigSeekComponent<Int, Real>
{
    // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
    IntVec highBoundaries; // Default values are 1 and b-1
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
    bool permax;   // When true, we retain only max-period generators. Default is false.
};


/**
 * Configuration for a seek from a .gen file.
 */
template<typename Int, typename Real> struct ConfigSeekGenFile : ConfigSeekComponent<Int, Real>
{
    string genFile;
};


/**
 * Factory implementation.
 */
template<typename Int, typename Real> ConfigSeekComponent<Int, Real>* ConfigSeekComponent<Int, Real>::create(GenType type)
{
    switch(type)
    {
        case MRG:
            return new ConfigSeekMRG<Int, Real>();

        case MWC:
            return new ConfigSeekMWC<Int, Real>();

        case GEN:
            return new ConfigSeekGenFile<Int, Real>();

        default:
            return nullptr;
    }
}


/**
 * Helper cast functions.
 */
template<typename Int, typename Real>
ConfigSeekMRG<Int, Real>*
asMRG(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekMRG<Int, Real>*>(ptr);
}

template<typename Int, typename Real>
ConfigSeekMWC<Int, Real>*
asMWC(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekMWC<Int, Real>*>(ptr);
}

template<typename Int, typename Real>
ConfigSeekGenFile<Int, Real>*
asGEN(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekGenFile<Int, Real>*>(ptr);
}


/**
 * This stores the information on the FOM that is used and how it is computed.
 */
template<typename Int, typename Real> struct ConfigMerit
{
    NTL::Vec<int64_t> t;
    LatticeTester::Weights* weights = nullptr;
    LatticeTester::Normalizer* norma = nullptr;
    ReducerBB<Int, Real> *red = 0;
    bool includeFirst = false;
    bool dualLattice = false;  

    // CW: Taken from old version -- probably not needed anymore 
    LatticeTester::NormaType normaType = LatticeTester::NONE;
    LatticeTester::CriterionType criterion = LatticeTester::SPECTRAL;
    LatticeTester::ReductionType reduction = LatticeTester::BB;
    LatticeTester::NormType norm = LatticeTester::L2NORM;

    // CW: Taken from old version -- maybe also not needed anymore
    #ifdef LATMRG_SEEK
        bool best = true;
        int num_gen = 0;
        double currentMerit = double(0);
        std::string construction = "RANDOM";
    #endif
    
    int max_gen = 10;

};


/**
 * This struct stores all the information required for a seek operation.
 */
template<typename Int, typename Real> struct ConfigSeek
{
    /* From the guide:
    * A) Types of generators to be searched:
    * Type => Variable genType
    * generators with C ≥ 1 components => Variable numComp
    * 
    * B) MRG component:
    * Type (MRG) => really necessary?
    * modulus m => Variable?
    * order k => Variable?
    * Rectangular region b=(b_1,...,b_k), c = (c_1,...,c_k) for multipliers 
    * Method: exhaustive or random
    * additional constraints: To be added later on
    * either search only for MRG components having maximal period, or ignore the period
    * 
    * C) MWC component
    * Type (MWC)
    * b
    * Rectangular region b=(b_1,...,b_k), c = (c_1,...,c_k) for a_0
    * maximal period
    * 
    * D) Matrix LCG component:
    * Not yet implemented
    * 
    * E) Reading generators from a file:
    * To be done later
    * 
    * F) Definition of output vectors and of the lattice Ls.    
    * analyze the lattice structure for vectors formed by groups of s successive values starting d values apart
    * 
    * G) Figures of merit
    * Set values necessary to define figure of Merit (currently only M is implemented)
    * 
    * H) Method of search
    * When examining a vector a, the program first checks if the maximal period conditions are satisfied, if this is required.
    * If a is not rejected by the maximal period test, then we move forward to the next MRG component and try all the vectors for that next 
    * component (by exhaustive or random search) and examine their combination with the currently examined multipliers for the previous components.
    * For each combined generator, we compute the FOM, but we interrupt the computation as soon as we find that this generator is not worth considering.
    * For this, the program always keeps lower and upper bounds on the FOM. These bounds are initialized at MinMerit and MaxMerit, and are updated whenever we can.
    * The current total execution (CPU) time is also checked before testing each new generator. When it exceeds the CPU time limit given in the data file, 
    * the search is aborted and the partial results are printed. 
    * One can provide a seed for the RNG that is used when the search is random. There is a default seed for when no seed is provided.
    * 
    * I) Output choices
    */

    GenType genType;   // Type of generator, currently MRG or MWC.
    long numComp;      // Number of components.  If > 1, we have a combined generator.
    int64_t maxdim;         // Maximal dimension of the lattice (new by CW: maybe not needed in the end)

    // List of configurations for the `numComp` components. They can be MRG or MWC.
    // The size of this list must be equal to numComp.
    // Each entry is a configuration for one component.
    vector<ConfigSeekComponent<Int, Real>*> genComponents; 

    bool leapBlocks;  // True iff we use lacunary indices with blocks at large leaps.
    Int leapSize;     // The leap size d when leapBlocks = true.
    Int blockSize;    // The block size s when leapBlocks = true.

    LatticeType latticeType; // The type of lattice that is constructed (see EnumTypes).
    // IntVect initState;  // The initial state, in case `latticeType = orbit`.
                           // The `orbit` case is currently not implemented. ***

    ConfigMerit<Int, Real> configFOM;  // The configuration for the FOM and its computation.

    SearchMethod searchMethod; // The method of search for the seek (exhaustive or random).
    long numRegions = 1;      // Number of regions for the search.
    vector<long> sizeRegions; // The sizes of the regions to be examined.
    long rngSeed = 123456;    // A seed for the RNG used in the search.
    long numRetained = 1;     // Number of generators that we want to retain.
    double timeLimit = 600;         // CPU time limit in seconds.

    OutputType outType;       // Type of output. 
    long verboseLevel = 1;    // Level of verbosity in the output.  Integer from 0 to 3.
    string resFile;           // If present, the results will be printed in file `resFile.res`.
    // string texFile;           // If present, the results will be printed in Latex in file `texFile.tex`.
    string genFile;           // If present, the results will be printed .gen format in file `genFile.gen`.
    bool showTimes = true;    // If true, prints CPU times for various operations.
    
   // Projection
   Projections* proj = nullptr;
   
   
    #ifdef LATMRG_SEEK
        bool progress = true; // Prints the program progress between generators
    #endif
    
    // Program variables - CW: probably not needed
    int detail = 1;


   // CW: Unclear what these variables are really used for    
    bool gen_set = false;
    bool test_set = false;
    bool proj_set = false;    

    
   // This is used to store the coefficients of a MRG and the information on how to
   // search for coefficients of a generator. This can have multiple different formats.
   std::vector<IntVec> coeff;


    // End inertion CW

    /**
     * Automatically creates the correct component
     * from conf.genType.
     */
    ConfigSeekComponent<Int, Real>* createComponent()
    {
        return ConfigSeekComponent<Int, Real>::create(genType);
    }

    ~ConfigSeek()
    {
        for(auto* ptr : genComponents)
            delete ptr;
    }
};

#endif

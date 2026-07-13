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

    /**
     * Factory function.
     * Creates the correct subclass automatically.
     */
    static ConfigSeekComponent<Int, Real>* create(GenType type);

    /**
     * Several virtual methods to get important information of the RNG.
     */
    virtual Int getModulus() const { return Int(0);}
    
    virtual long getOrder() const { return 0;}

    virtual Int getNoMultipliers() const { return Int(0);}
    
    virtual Int getLowBoundary(int i) const { return Int(0);}

    virtual Int getHighBoundary(int i) const { return Int(0);}
    
    virtual int64_t getRandomBits(int i) const { return 0;}
    
    virtual long getPowMod() const { return long(0);}

    virtual bool onlyMaxPeriod() const { return false;} 

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
    bool permaxPrime = false;  // When true, m must be prime and we retain only max-period MRGs (default = false)
    DecompType howFactorh; // Indicates how to factorize h=(m-1)/2
    IntVec factorsh;      // The factors of h, repeated according to multiplicity, when known
    DecompType howFactor;  // Indicates how to factorize r=(m^k-1)/(m-1)
    IntVec factorsr;      // The factors of r, repeated according to multiplicity, when known
    bool permaxPow2;   // True iff m is a power of 2, b=2, r=0, k=1, and we want max period for the LCG
                       // Default is false

    Int getModulus() const override { return modulus; }
    
    long getOrder() const override { return order; }

    Int getNoMultipliers() const {
      Int total(1);
      for (long i = 1; i <= order; i++) {
          total *= (highBoundaries[i] - lowBoundaries[i] + 1);
      }
      return total;
    }
    
    Int getLowBoundary(int i) const { return lowBoundaries(i);}
    
    Int getHighBoundary(int i) const { return highBoundaries(i);}

    bool onlyMaxPeriod() const { return permaxPrime;} 
};


/**
 * Configuration for a MWC component.
 */
template<typename Int, typename Real> struct ConfigSeekMWC : ConfigSeekComponent<Int, Real>
{
    // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
    IntVec highBoundaries; // Default values are 1 and b-1
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<int64_t> randomBits; // For random searches, the nonzero coefficients `a_j` are generated randomly with `e_j` random bits
    int64_t numaj; // number of positive coefficients 
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
    bool permax = false;   // When true, we retain only max-period generators. Default is false.

    
    long getOrder() const override { return order; }

    Int getNoMultipliers() const {
      Int total(1);
      for (long i = 1; i <= order; i++) {
          total *= (highBoundaries[i] - lowBoundaries[i] + 1);
      }
      return total;
    }
    
    Int getLowBoundary(int i) const { return lowBoundaries(i);}
    
    Int getHighBoundary(int i) const { return highBoundaries(i);}
    
    int64_t getRandomBits(int i) const { return randomBits[i];}

    long getPowMod() const { return powMod; }
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

        default:
            return nullptr;
    }
}


/**
 * Helper cast functions.
 */
template<typename Int, typename Real> ConfigSeekMRG<Int, Real>* asMRG(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekMRG<Int, Real>*>(ptr);
}

template<typename Int, typename Real> ConfigSeekMWC<Int, Real>* asMWC(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekMWC<Int, Real>*>(ptr);
}

template<typename Int, typename Real> ConfigSeekGenFile<Int, Real>* asGEN(ConfigSeekComponent<Int, Real>* ptr)
{
    return dynamic_cast<ConfigSeekGenFile<Int, Real>*>(ptr);
}


/**
 * ConfigMerit stores the information on the FOM that is used and how it is computed.
 */
template<typename Int, typename Real> struct ConfigMerit
{
    NTL::Vec<int64_t> t;
    LatticeTester::Weights* weights = nullptr;
    LatticeTester::Normalizer* norma = nullptr;
    ReducerBB<Int, Real> *red = nullptr;
    bool includeFirst = true;
    bool dualLattice = false; 
    LatticeTester::NormType norm = LatticeTester::L2NORM;
    
    #ifdef LATMRG_SEEK
        bool best = true;
        int num_gen = 0;
        double currentMerit = double(0);
    #endif
    
    int max_gen = 10; // Maximal number of generators to be tested

};


/**
 * This struct stores all the information required for a general seek operation.
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
    * Rectangular region b=(b_1,...,b_k), c = (c_1,...,c_k) for multipliers (*)
    * Method: exhaustive or random (only exhaustive so far)
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
    int64_t maxdim;    // Maximal dimension of the lattice

    // List of configurations for the `numComp` components. They can be MRG or MWC.
    // The size of this list must be equal to numComp.
    // Each entry is a configuration for one component.
    vector<ConfigSeekComponent<Int, Real>*> genComponents; 

    LatticeType latticeType; // The type of lattice that is constructed (see EnumTypes).
    
    ConfigMerit<Int, Real> configFOM;  // The configuration for the FOM and its computation.

    SearchMethod searchMethod; // The method of search for the seek (exhaustive or random).
    double timeLimit = 600;         // CPU time limit in seconds.

    #ifdef LATMRG_SEEK
        bool progress = true; // Prints the program progress between generators
    #endif

    /**
     * Automatically creates the correct component
     * from conf.genType.
     */
    void createComponent()
    {
        genComponents.resize(numComp);
        for (long i = 0; i < numComp; i++) {
            genComponents[i] = ConfigSeekComponent<Int, Real>::create(genType);
        }
    }

    ~ConfigSeek()
    {
        //for(auto* ptr : genComponents)
        //    delete ptr;
    }
};

#endif

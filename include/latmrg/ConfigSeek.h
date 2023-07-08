#ifndef LATMRG_CONFIGSEEK_H
#define LATMRG_CONFIGSEEK_H

#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/Reducer.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/ComboLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/Projections.h"
#include "latmrg/FigureOfMerit.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/Projections.h"
#include "latmrg/Test.h"

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace std;
using namespace LatMRG;
using namespace LatticeTester::Random;

typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

/**
 * This `struct` stores all the information required for a seek operation.
 * The information is organized in the same way as in a XML file for a seek.
 * The variables names also correspond to the tag names.
 */
template<typename Int, typename Real> struct ConfigSeek {

    GenType genType;   // Type of generator, currently MRG or MWC.
    long numComp;      // Number of components.  If > 1, we have a combined generator.

    // List of configurations for the `numComp` components. They can be MRG or MWC.
    // The size of this list must be equal to numComp.
    // Each entry is a configuration for one component.
    vector<ConfigSeekComponent<Int>*> genComponents;

    bool leapBlocks;  // True iff we use lacunary indices with blocks at large leaps.
    Int leapSize;     // The leap size d when leapBlocks = true.
    Int blockSize;    // The block size s when leapBlocks = true.

    LatticeType latticeType; // The type of lattice that is constructed (see EnumTypes).
    // IntVect initState;  // The initial state, in case `latticeType = orbit`.
                           // The `orbit` case is currently not implemented. ***

    ConfigMerit configFOM;  // The configuration for the FOM and its computation.

    SearchMethod searchMethod; // The method of search for the seek (exhaustive or random).
    long numRegions = 1;      // Number of regions for the search.
    vector<long> sizeRegions; // The sizes of the regions to be examined.
    long rngSeed = 123456;    // A seed for the RNG used in the search.
    long numRetained = 1;     // Number of generators that we want to retain.
    double timeLimit;         // CPU time limit in seconds.

    OutputType outType;       // Type of output.
    long verboseLevel = 1;    // Level of verbosity in the output.  Integer from 0 to 3.
    string resFile;           // If present, the results will be printed in file `resFile.res`.
    // string texFile;           // If present, the results will be printed in Latex in file `texFile.tex`.
    string genFile;           // If present, the results will be printed .gen format in file `genFile.gen`.
    bool showTimes = true;    // If true, prints CPU times for various operations.
};



/**
 * Abstract structure to store the information required for a component in a seek.
 */
template<typename Int, typename Real> abstract struct ConfigSeekComponent;


/**
 * This `struct` stores the information required for a seek for a MRG component.
 */
template<typename Int, typename Real> struct ConfigSeekMRG : ConfigSeekComponent {
    Int modulus;   // The modulus m = b^e + r
    Int base;      // The base b
    Int exponent;  // The exponent e
    Int rest;      // The rest r
    long order;    // The order k
    Intvect lowBoundaries;  // The low and high boundaries for the coefficients a_j
    Intvect highBoundaries; // Default values are 1 and m-1
    bool appFact;   // true iff we enforce the approximate factoring condition (default = false)
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
    bool permaxPrime;  // When true, m must be prime and we retain only max-period MRGs (default = false)
    DecompType howFactorh; // Indicates how to factorize h=(m-1)/2
    IntVect factorsh;      // The factors of h, repeated according to multiplicity, when known
    DecompType howFactor;  // Indicates how to factorize r=(m^k-1)/(m-1)
    IntVect factorsr;      // The factors of r, repeated according to multiplicity, when known
    bool permaxPow2;   // True iff m is a power of 2, b=2, r=0, k=1, and we want max period for the LCG
                       // Default is false
};

/**
 * This stores the information required for a seek for a MWC component.
 */
template<typename Int, typename Real> struct ConfigSeekMWC : ConfigSeekComponent {
    // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    Intvect lowBoundaries;  // The low and high boundaries for the coefficients a_j
    Intvect highBoundaries; // Default values are 1 and b-1
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
    bool permax;   // When true, we retain only max-period generators. Default is false.
};


/**
 * This stores the information required for a seek from a .gen file.
 */
template<typename Int, typename Real> struct ConfigSeekGenFile : ConfigSeekComponent {
    string genFile;
};


/**
 * This stores the information on the FOM that is used and how it is computed.
 */
template<typename Int, typename Real> struct ConfigMerit {
    MeritType meritType;  // The type of FOM that is used.
    long meritd;   // value of d in the FOM.
    vector<long> meritDims;
    bool maximize = true;  // true iff we maximize the FOM.
    double lowbound = 0.0;   // Lower bound on the FOM.  Default is 0.
    double highbound = 1.0;   // Upper bound on the FOM.  Default is 1.
    bool fomInDual = true; // True iff the FOM is based on the dual lattice.
    NormType norm = L2Norm;  // The norm used to measure the vector lengths, Default is `L2Norm`.
    NormaType normalizer = BestLat;  // Type of normalization that is used.

    PreReductionType prered;  // Type of pre-reduction of the basis.
    double delta = 0.999999;  // For LLL and BKZ.
    blockSizeBKZ = 10;        // Block size for BKZ.
	BBDecompType decompBB;    // Type of decomposition for the BB algorithm (Cholesky or triangular)
	ProjConstructType projConstruct;  // Method to construct the projections (LLL or triangular).
    long maxnodesBB = 1 << 8; // Max number of nodes to be examined by the BB.
    };




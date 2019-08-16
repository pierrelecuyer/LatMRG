typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

// For period tests
DecompType decompm1 = DECOMP, decompr = DECOMP;
std::string filem1, filer;

// Data file read parameters
GenType type = MRG;
// Type of figure of merit
LatticeTester::NormaType normaType = LatticeTester::NONE;
LatticeTester::CriterionType criterion = LatticeTester::SPECTRAL;
LatticeTester::PreReductionType reduction = LatticeTester::FULL;
bool use_dual = true;
#ifdef LATMRG_SEEK
bool best = true;
int num_gen = 0;
Dbl currentMerit = Dbl(0);
std::string construction = "RANDOM";
bool gen_set = false;
bool test_set = false;
bool proj_set = false;
#endif
// Projection
Projections* proj;

double timeLimit = 600;
int max_gen = 10;

// MRG Specific parameters
IntVec mult; // MRG multipliers
IntVec coeff;

// MWC Specific parameters
Int b; // modulo of MWC recurence

// MMRG Specific parameters
IntMat matrix;

// Combo specific parameters
#ifdef LATMRG_SEEK
int num_comp;
std::vector<std::int64_t> comb_order; // k for MRG and MMRG
IntVec comb_modulo; // m for MRG and MMRG
// modulo is basis^exponent+rest
std::vector<std::int64_t> comb_basis;
std::vector<std::int64_t> comb_exponent;
std::vector<std::int64_t> comb_rest;
std::vector<bool> comb_period;
std::vector<MRGComponent<Int>*> comb_fact;
#endif

// Shared components names
std::int64_t order; // k for MRG and MMRG
Int modulo; // m for MRG and MMRG
// modulo is basis^exponent+rest
std::int64_t basis;
std::int64_t exponent;
std::int64_t rest;
bool period = true;

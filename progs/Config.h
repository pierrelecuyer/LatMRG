typedef Integ Int;
typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

// Type of figure of merit
LatticeTester::NormaType normaType = LatticeTester::NONE;
LatticeTester::CriterionType criterion = LatticeTester::SPECTRAL;
LatticeTester::PreReductionType reduction = LatticeTester::FULL;
LatticeTester::NormType norm = LatticeTester::L2NORM;

bool use_dual = true;
#ifdef LATMRG_SEEK
bool best = true;
int num_gen = 0;
Dbl currentMerit = Dbl(0);
std::string construction = "RANDOM";
#endif
// Projection
Projections* proj;

int max_gen = 10;

// MWC Specific parameters
Int b; // modulo of MWC recurence

// MRGComponent stores the stuff we might want to know, such as the modulo
// the order and even the coefficients
int num_comp = 1;
std::vector<string> search_mode;
std::vector<bool> period;
std::vector<MRGComponent<Int>*> fact;
// This is used to store the coefficients of a MRG and the information on how to
// search for coefficients of a generator. This can have multiple different formats.
std::vector<IntVec> coeff;

bool gen_set = false;
bool test_set = false;
bool proj_set = false;

// Program variables
double timeLimit = 600;
int detail = 1;
#ifdef LATMRG_SEEK
bool progress = true; // Prints the program progress between generators
#endif

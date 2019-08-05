#include <cstring>
#include <string>
#include <iostream>
#include <fstream>

#include "tinyxml2.h"

//#include "TestLattice.h"

#define MRGLATTICE_MAIN_EXEC
#include "ExecCommon.h"

#include "Seek.h"

using namespace LatMRG;

bool options = true;
std::string exec_mode;

// Output stuff
std::ofstream fout;
std::ostream* out(&std::cout);

// Prints the program usage
void print_help() {
  *out << "Usage: MRGLattice [OPTIONS] MODE PARAMETERS\n\n";
  *out << "Available options:\n\t" << "-h, --help\n\t\tPrint this"
   " message.\n\t-o, --output=FILE\n\t\tOutput to FILE\n\t"
   "-t, --types=[LD,ZD,ZR]\n\t\tUse corresponding types. (default: ZD)\n\n";
  *out << "Available modes:\n\t" << "lat\n\t\t" << "Study the"
    " lattice of a MRG generator.\n\t" << "seek\n\t\t" << "Search for"
    " new MRG generators.\n\t" << "mk\n\t\t" << "Find compatible "
    "modulus and order for MRG generators.\n\t" << "period\n\t\t" <<
    "Test the period of a MRG generator.\n\n";
  *out << "lat parameters:\n\t" << "    --gentype=[MRG,MWC,COMBO,MMRG]\n\t\t"
    "Use generator type GenType. (default: MRG)\n\t"
    "    --criterion=[LENGTH,SPECTRAL, BEYER]\n\t\t"
    "Criterion for the figure of merit. (default: SPECTRAL)\n\t"
    "    --dual=[true,false]\n\t\tCompute on the dual or no. (default: true)\n\t"
    "    --reduction\n\t" << "    --normalizer\n\t" << "    --time\n\t"
    "    --projections\n\t" << "-a, --vector\n\t" << "-m, --modulo\n\t"
    "-k, --order\n";
}

//==============================================================================
//===== Miscellanous functions to move somewhere else
//==============================================================================

int toVectString(const char* str, IntVec& vect, int length) {
  int j = 0, old = 0;;
  for (int i = 0; i < length; i++) {
    while (str[j] != ' ' && str[j] != '\0') j++;
    if (str[j] == '\0' && j <= old) {
      std::cerr << "Wrong number of elements in attribute 'random' of tag 'method'.\n";
      return 1;
    }
    if (j <= old) {
      i--;
      continue;
    }
    char to[j-old];
    strncpy(to, str+old, j-old);
    old = j = j+1;
    NTL::conv(vect[i], to);
  }
  return 0;
}

//==============================================================================
//===== Reading tags used in multiple configureations
//==============================================================================

// Read all info of a MRG prints error message if essential info is missing
int readMRG(tinyxml2::XMLNode* current) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("modulo");
  if (node) {
    NTL::conv(p, node->Attribute("p"));
    NTL::conv(r, node->Attribute("r"));
    NTL::conv(e, node->Attribute("e"));
    m = NTL::power(p, e) + r;
  } else {
    std::cerr << "No 'modulo' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(k, node->Attribute("k"));
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("method");
  if (node) {
    auto value = node->FirstAttribute();
    if (value) {
      if (!strcmp(value->Name(), "random")) {
        coeff.SetLength(k);
        if (toVectString(value->Value(), coeff, k)) return 1;
      } else if (!strcmp(value->Name(), "pow2")) {
        coeff.SetLength(2*k);
      } else {
        std::cerr << "Invalid attribute in tag 'method'.\n";
        return 1;
      }
    } else {
      std::cerr << "Tag 'method' has no attribute.\n";
      return 1;
    }
  } else {
    node = current->FirstChildElement("coefficients");
    if (node) {
      // do stuff
    } else {
      std::cerr << "No way to set coefficients in 'mrg' tag. Add 'method' or 'coefficients' tag.\n";
      return 1;
    }
  }

  return 0;
}

//==============================================================================

// Reads parameters for a spectral test. Will print errors but not halt the
// program since all values have a default
int readSpectral(tinyxml2::XMLNode* current) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("reduction");
  if (node) {
    if (toRedString(red_type, node->Attribute("method"))) std::cerr << "Invalid"
      "/inexistant 'method' attribute in 'reduction' tag.\n";
  }

  node = current->FirstChildElement("norma");
  if (node) {
    if (toNormaString(norma_type, node->Attribute("lizer"))) std::cerr << "Invalid"
      "/inexistant 'lizer' attribute in 'norma' tag.\n";
  }
  return 0;
}

//==============================================================================
//===== Dispatching reused tags to the right functions
//==============================================================================

int tryGen(tinyxml2::XMLNode* current) {
  //tinyxml2::XMLNode* node;
  if (!strcmp(current->Value(), "mrg")) {
    gen_type = MRG;
    if (readMRG(current)) {
      std::cerr << "'mrg' tag is not properly parametrized.\n";
      return 1;
    }
    gen_set = true;
  } else if (!strcmp(current->Value(), "mmrg")) {
    std::cout << current->Value() << "\n";
  } else if (!strcmp(current->Value(), "mwc")) {
    std::cout << current->Value() << "\n";
  } else if (!strcmp(current->Value(), "combo")) {
    std::cout << current->Value() << "\n";
  } else return 0;
  return 1;
}

//==============================================================================

int tryTest(tinyxml2::XMLNode* current) {
  //tinyxml2::XMLNode* node;
  if (!strcmp(current->Value(), "spectral")) {
    crit_type = LatticeTester::SPECTRAL;
    readSpectral(current);
  } else if (!strcmp(current->Value(), "beyer")) {
    std::cout << current->Value() << "\n";
  } else if (!strcmp(current->Value(), "shortest")) {
    std::cout << current->Value() << "\n";
  } else return 0;
  return 1;
}

//==============================================================================
//===== Different uses reding functions
//==============================================================================

void readMK(tinyxml2::XMLNode* current) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!strcmp(node->Value(), "bounds")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "k")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "safe")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "factor")) {
      std::cout << node->Value() << "\n";
    }
    node = node->NextSibling();
  }
}

//==============================================================================

void readPeriod(tinyxml2::XMLNode* current) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!strcmp(node->Value(), "bounds")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "k")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "safe")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "factor")) {
      std::cout << node->Value() << "\n";
    }
    node = node->NextSibling();
  }
}

//==============================================================================

void readSeek(tinyxml2::XMLNode* current) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (tryGen(node)) {}
    else if (tryTest(node)) {}
    else if (!strcmp(node->Value(), "bounds")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "k")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "safe")) {
      std::cout << node->Value() << "\n";
    } else if (!strcmp(node->Value(), "factor")) {
      std::cout << node->Value() << "\n";
    }
    node = node->NextSibling();
  }
  conf.type = gen_type;
  conf.normaType = norma_type;
  conf.criterion = crit_type;
  conf.reduction = red_type;
  conf.best = best;
  conf.timeLimit = time_limit;
  conf.max_gen = max_gen;
  conf.order = k;
  conf.basis = p;
  conf.modulo = m;
  conf.exponent = e;
  conf.rest = r;
  conf.period = period;
}

//==============================================================================
//===== Main reading function
//==============================================================================

void readFile(const char* filename) {
  tinyxml2::XMLDocument doc;
  doc.LoadFile(filename);
  tinyxml2::XMLNode* current;
  current = doc.FirstChild();
  if (!strcmp(current->Value(), "mk")) {
    readMK(current);
    // fill mk
  } else if (!strcmp(current->Value(), "period")) {
    readPeriod(current);
    // fill period
  } else if (!strcmp(current->Value(), "seek")) {
    readSeek(current);
    Seek();
    // fill seek
  } else if (!strcmp(current->Value(), "lattest")) {
    std::cout << 1 << "\n";
    // fill lattest
  }
}

//==============================================================================

int main(int argc, char** argv) {
  if (argc < 2) {
    print_help();
    return 1;
  }
  for (int i = 1; i < argc; i++) {
    readFile(argv[i]);
  }
  //if (exec_mode == "lat") {
  //  return TestLattice(gen_type, crit_type, dual, red_type, norma_type,
  //      time_limit, detail, period, *proj, gen_string);
  //} else if (exec_mode == "seek") {
  //} else if (exec_mode == "mk") {
  //} else if (exec_mode == "period") {
  //}
  return 0;
}

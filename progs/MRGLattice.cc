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

/**
 * Fills the first `length` values of `vect` with elements in `str`.
 * This functions assumes `vect` is at least `length` elements long, that its
 * elements can be accessed with `operator[]` and that there exists an overload
 * `NTL::conv(vect_type, char*)` between their type and a `C` string.
 * */
template<typename Vec>
int toVectString(const char* str, Vec& vect, int length) {
  int j = 0, old = 0;
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
template<typename Conf>
int readMRG(tinyxml2::XMLNode* current, Conf& conf) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("modulo");
  if (node) {
    NTL::conv(conf.basis, node->Attribute("basis"));
    NTL::conv(conf.rest, node->Attribute("rest"));
    NTL::conv(conf.exponent, node->Attribute("exponent"));
    conf.modulo = NTL::power(conf.basis, conf.exponent) + conf.rest;
  } else {
    std::cerr << "No 'modulo' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(conf.order, node->Attribute("k"));
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("method");
  if (node) {
    auto value = node->FirstAttribute();
    if (value) {
      if (!strcmp(value->Name(), "random")) {
        conf.coeff.SetLength(conf.order);
        if (toVectString(value->Value(), conf.coeff, conf.order)) return 1;
      } else if (!strcmp(value->Name(), "pow2")) {
        conf.coeff.SetLength(2*conf.order);
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

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf))
    std::cerr << "Non critical error in 'period' tag of 'mrg' tag.\n";

  return 0;
}

//==============================================================================

// Read all info of a MMRG prints error message if essential info is missing
template<typename Conf>
int readMMRG(tinyxml2::XMLNode* current, Conf& conf) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("modulo");
  if (node) {
    NTL::conv(conf.basis, node->Attribute("basis"));
    NTL::conv(conf.rest, node->Attribute("rest"));
    NTL::conv(conf.exponent, node->Attribute("exponent"));
    conf.modulo = NTL::power(conf.basis, conf.exponent) + conf.rest;
  } else {
    std::cerr << "No 'modulo' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(conf.order, node->Attribute("k"));
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf))
    std::cerr << "Non critical error in 'period' tag of 'mmrg' tag.\n";

  return 0;
}

//==============================================================================

// Reads parameters for a spectral test. Will print errors but not halt the
// program since all values have a default
template<typename Conf>
int readSpectral(tinyxml2::XMLNode* current, Conf& conf) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("reduction");
  if (node) {
    if (toRedString(conf.reduction, node->FirstAttribute()->Value())) std::cerr << "Invalid"
      "/inexistant attribute in 'reduction' tag.\n";
  }

  node = current->FirstChildElement("norma");
  if (node) {
    if (toNormaString(conf.normaType, node->FirstAttribute()->Value())) std::cerr << "Invalid"
      "/inexistant attribute in 'norma' tag.\n";
  }

  node = current->FirstChildElement("dual");
  if (node) {
    conf.use_dual = node->FirstAttribute()->BoolValue();
  }
  return 0;
}

//==============================================================================

template<typename Conf>
int readProj(tinyxml2::XMLNode* current, Conf& conf) {
  int minDim, numProj;
  std::vector<int> projDim;
  auto node = current->FirstChildElement("min");
  if (!node || !(node->FirstAttribute())) {
    std::cerr << "Missing/incorrect child node 'min' in 'proj' tag.\n";
    return 1;
  }
  minDim = node->FirstAttribute()->IntValue();
  node = current->FirstChildElement("num");
  if (!node || !(node->FirstAttribute())) {
    std::cerr << "Missing/incorrect child node 'num' in 'proj' tag.\n";
    return 1;
  }
  numProj = node->FirstAttribute()->IntValue();
  projDim.resize(numProj);
  node = current->FirstChildElement("dim");
  if (!node || !(node->FirstAttribute())) {
    std::cerr << "Missing/incorrect child node 'dim' in 'proj' tag.\n";
    return 1;
  }
  if(!toVectString(node->FirstAttribute()->Value(), projDim, numProj)) {
    for (int i = 0; i < numProj; i++){
      if (conf.max_dim < projDim[i]) conf.max_dim = projDim[i];
    }
    conf.proj = new Projections(numProj, minDim, projDim);
    conf.proj_set = true;
    return 0;
  }
  return 1;
}

//==============================================================================

// Reads info on how to check for generator maximal period.
template<typename Conf>
int readGenPer(tinyxml2::XMLNode* current, Conf& conf) {
  conf.period = current->ToElement()->BoolAttribute("check");
  int err = 0;
  if (conf.period) {
    auto node = current->FirstChildElement("m1");
    if (node) {
      auto method = node->Attribute("method");
      if (method && toDecomString(conf.decompm1, method)) {
        std::cerr << "Invalid 'method' attribute value in 'm1' tag.\n";
        err = 1;
      }
      auto name = node->Attribute("file");
      if (name) conf.filem1 = name;
    }
    node = current->FirstChildElement("r");
    if (node) {
      auto method = node->Attribute("method");
      if (method && toDecomString(conf.decompr, method)) {
        std::cerr << "Invalid 'method' attribute value in 'r' tag.\n";
        err = 1;
      }
      auto name = node->Attribute("file");
      if (name) conf.filer = name;
    }
  }
  return err;
}

//==============================================================================
//===== Dispatching reused tags to the right functions
//==============================================================================

template<typename Conf>
int tryGen(tinyxml2::XMLNode* current, Conf& conf) {
  if (!strcmp(current->Value(), "mrg")) {
    conf.type = MRG;
    if (readMRG(current, conf)) {
      std::cerr << "'mrg' tag is not properly parametrized.\n";
      return 1;
    }
    conf.gen_set = true;
  } else if (!strcmp(current->Value(), "mmrg")) {
    conf.type = MMRG;
    if (readMMRG(current, conf)) {
      std::cerr << "'mmrg' tag is not properly parametrized.\n";
      return 1;
    }
    conf.gen_set = true;
  } else if (!strcmp(current->Value(), "mwc")) {
    std::cout << current->Value() << "\n";
  } else if (!strcmp(current->Value(), "combo")) {
    std::cout << current->Value() << "\n";
  } else return 0;
  return 1;
}

//==============================================================================

template<typename Conf>
int tryTest(tinyxml2::XMLNode* current, Conf& conf) {
  //tinyxml2::XMLNode* node;
  if (!strcmp(current->Value(), "spectral")) {
    conf.criterion = LatticeTester::SPECTRAL;
    readSpectral(current, conf);
    conf.test_set = true;
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

template<typename Conf>
void readSeek(tinyxml2::XMLNode* current, Conf& conf) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!conf.gen_set && tryGen(node, conf)) {
    } else if (!conf.test_set && tryTest(node, conf)) {
    } else if (!conf.proj_set && !strcmp(node->Value(), "proj")) {
      readProj(node, conf);
    } else if (!strcmp(node->Value(), "time")) {
      conf.timeLimit = node->ToElement()->FirstAttribute()->DoubleValue();
    } else if (!strcmp(node->Value(), "num_gen")) {
      conf.max_gen = node->ToElement()->FirstAttribute()->IntValue();
    } else if (!strcmp(node->Value(), "best")) {
      conf.best = node->ToElement()->FirstAttribute()->BoolValue();
      conf.currentMerit = conf.best?0.0:1.0;
    }
    node = node->NextSibling();
  }
}

//==============================================================================

template<typename Conf>
void readTest(tinyxml2::XMLNode* current, Conf& conf) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!conf.gen_set && tryGen(node, conf)) {
      std::cout << 1 << "\n";
    } else if (!conf.test_set && tryTest(node, conf)) {
      std::cout << 2 << "\n";
    } else if (!conf.proj_set && !strcmp(node->Value(), "proj")) {
      std::cout << 3 << "\n";
      readProj(node, conf);
    } else if (!strcmp(node->Value(), "time")) {
      std::cout << 4 << "\n";
      conf.timeLimit = node->ToElement()->FirstAttribute()->DoubleValue();
    } else if (!strcmp(node->Value(), "num_gen")) {
      std::cout << 5 << "\n";
      conf.max_gen = node->ToElement()->FirstAttribute()->IntValue();
    }
    node = node->NextSibling();
  }
}


//==============================================================================
//===== Main reading function
//==============================================================================

int readFile(const char* filename) {
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
    if (types == "LD") {
      SeekMain<std::int64_t, double> prog;
      readSeek(current, prog.conf);
      return prog.Seek();
    } else if (types == "ZD") {
      SeekMain<NTL::ZZ, double> prog;
      readSeek(current, prog.conf);
      return prog.Seek();
    } else if (types == "ZR") {
      SeekMain<NTL::ZZ, NTL::RR> prog;
      readSeek(current, prog.conf);
      return prog.Seek();
    }
  } else if (!strcmp(current->Value(), "lattest")) {
    //LatTest<NTL::ZZ, NTL::RR> prog;
    //readTest(current, prog.conf);
    //return prog.TestLat();
  }
  return 0;
}

//==============================================================================

int main(int argc, char** argv) {
  if (argc < 2) {
    print_help();
    return 1;
  }
  for (int i = 1; i < argc; i++) {
    if (readFile(argv[i])) {
      std::cerr << "File " << argv[i] << " exited with errors.\n";
    }
  }
  return 0;
}

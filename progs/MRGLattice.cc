#include <cstring>
#include <string>
#include <iostream>
#include <fstream>

#include "tinyxml2.h"

//#include "TestLattice.h"

#define MRGLATTICE_MAIN_EXEC
#include "ExecCommon.h"

#include "Seek.h"
#include "LatTest.h"
#include "FindMK.h"

using namespace LatMRG;

bool options = true;
std::string exec_mode;

// Output stuff
std::ofstream fout;
std::ostream* out(&std::cout);

// Prints the program usage
void print_help() {
  *out << "Usage: MRGLattice file1 file2 ...\n";
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
    char to[j-old+1];
    strncpy(to, str+old, j-old);
    to[j-old] = '\0';
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
int readMRG(tinyxml2::XMLNode* current, Conf& conf, int i) {
  tinyxml2::XMLElement* node;
  std::cout << 1 << "\n";

  int exponent, order;
  typename Conf::Int /*modulo,*/ basis, rest;

  std::string filem1 = "./tempm1" + i + std::to_string(rand());
  std::string filer = "./tempr" + i + std::to_string(rand());
  DecompType decompm1 = DECOMP, decompr = DECOMP;

  std::cout << 2 << "\n";
  node = current->FirstChildElement("modulo");
  if (node) {
    if (!node->Attribute("basis")) {
      std::cerr << "No `basis` attribute in modulo tag.\n";
      return 1;
    }
    NTL::conv(basis, node->Attribute("basis"));
    if (node->Attribute("exponent")) NTL::conv(exponent, node->Attribute("exponent"));
    else exponent = 1;
    if (node->Attribute("rest")) NTL::conv(rest, node->Attribute("rest"));
    else rest = 0;
    //modulo = NTL::power(basis, exponent) + rest;
  } else {
    std::cerr << "No 'modulo' tag in 'mrg' tag.\n";
    return 1;
  }
  std::cout << 3 << "\n";

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(order, node->FirstAttribute()->IntValue());
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }
  std::cout << 4 << "\n";

  node = current->FirstChildElement("method");
  if (node) {
    auto value = node->FirstAttribute();
    if (value) {
      if (!strcmp(value->Name(), "random")) {
        conf.coeff[i].SetLength(order);
        if (toVectString(value->Value(), conf.coeff[i], order)) return 1;
      } else if (!strcmp(value->Name(), "pow2")) {
        conf.coeff[i].SetLength(2*order);
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
      conf.coeff[i].SetLength(order);
      toVectString(node->FirstAttribute()->Value(), conf.coeff[i], order);
    } else {
      std::cerr << "No way to set coefficients in 'mrg' tag. Add 'method' or 'coefficients' tag.\n";
      return 1;
    }
  }
  std::cout << 5 << "\n";

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf, i, filem1, filer, decompm1, decompr))
    std::cerr << "Non critical error in 'period' tag of 'mrg' tag.\n";
  std::cout << 6 << "\n";

  auto comp = new MRGComponent<typename Conf::Int>(basis, exponent, rest, order,
      decompm1, filem1.c_str(), decompr, filer.c_str());
  comp->set_type(MRG);

  conf.fact.push_back(comp);

  return 0;
}

//==============================================================================

// Read all info of a MMRG prints error message if essential info is missing
template<typename Conf>
int readMMRG(tinyxml2::XMLNode* current, Conf& conf, int i) {
  tinyxml2::XMLElement* node;

  int exponent, order;
  typename Conf::Int /*modulo,*/ basis, rest;

  std::string filem1 = "./tempm1" + i + std::to_string(rand());
  std::string filer = "./tempr" + i + std::to_string(rand());
  DecompType decompm1 = DECOMP, decompr = DECOMP;

  node = current->FirstChildElement("modulo");
  if (node) {
    NTL::conv(basis, node->Attribute("basis"));
    NTL::conv(rest, node->Attribute("rest"));
    NTL::conv(exponent, node->Attribute("exponent"));
    //conf.modulo[i] = NTL::power(conf.basis[i], conf.exponent[i]) + conf.rest[i];
  } else {
    std::cerr << "No 'modulo' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(order, node->FirstAttribute()->Value());
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf, i, filem1, filer, decompm1, decompr))
    std::cerr << "Non critical error in 'period' tag of 'mmrg' tag.\n";

  auto comp = new MRGComponent<typename Conf::Int>(basis, exponent, rest, order,
      decompm1, filem1.c_str(), decompr, filer.c_str());
  comp->set_type(MMRG);

  conf.fact.push_back(comp);

  return 0;
}

//==============================================================================

template<typename Conf>
int readMWC(tinyxml2::XMLNode* current, Conf& conf) {
  return 1;
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

// Reads parameters for shortest vector length computation
template<typename Conf>
int readLength(tinyxml2::XMLNode* current, Conf& conf) {
  tinyxml2::XMLElement* node;
  node = current->FirstChildElement("reduction");
  if (node) {
    if (toRedString(conf.reduction, node->FirstAttribute()->Value())) std::cerr << "Invalid"
      "/inexistant attribute in 'reduction' tag.\n";
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
    return 0;
  }
  return 1;
}

//==============================================================================

// Reads info on how to check for generator maximal period.
template<typename Conf>
int readGenPer(tinyxml2::XMLNode* current, Conf& conf, int i,
    std::string& filem1, std::string& filer, DecompType& decompm1, DecompType& decompr) {
  conf.period[i] = current->ToElement()->FirstAttribute()->BoolValue();
  int err = 0;
  if (conf.period[i]) {
    auto node = current->FirstChildElement("m1");
    if (node) {
      auto method = node->Attribute("method");
      if (method && toDecomString(decompm1, method)) {
        std::cerr << "Invalid 'method' attribute value in 'm1' tag.\n";
        err = 1;
      }
      auto name = node->Attribute("file");
      if (name) filem1 = name;
    }
    node = current->FirstChildElement("r");
    if (node) {
      auto method = node->Attribute("method");
      if (method && toDecomString(decompr, method)) {
        std::cerr << "Invalid 'method' attribute value in 'r' tag.\n";
        err = 1;
      }
      auto name = node->Attribute("file");
      if (name) filer = name;
    }
  }
  if (err) return 1;
  return 0;
}

//==============================================================================

template<typename Conf>
int readGenerator(tinyxml2::XMLNode* current, Conf& conf) {
  std::cout << 0 << "\n";
  auto tmp_node = current->FirstChildElement("numcomp");
  if (tmp_node) {
    conf.num_comp = tmp_node->FirstAttribute()->IntValue();
  }

  conf.coeff.resize(conf.num_comp);
  conf.period.resize(conf.num_comp);

  srand(time(NULL));

  std::cout << 1 << "\n";
  auto node = current->FirstChild();
  int count = 0;
  while (node && (count < conf.num_comp)) {
  std::cout << 2 << "\n";
    if (!strcmp(node->Value(), "mrg")) {
  std::cout << 3 << "\n";
      if (readMRG(node, conf, count)) {
        std::cerr << "'mrg' tag is not properly parametrized.\n";
        return 1;
      }
      count++;
    } else if (!strcmp(node->Value(), "mmrg")) {
      if (readMMRG(node, conf, count)) {
        std::cerr << "'mmrg' tag is not properly parametrized.\n";
        return 1;
      }
      count++;
    } else if (!strcmp(node->Value(), "mwc")) {
      std::cout << node->Value() << "\n";
      count++;
    }
    node->NextSibling();
  std::cout << 4 << "\n";
  }
  std::cout << 5 << "\n";

  if (count != conf.num_comp || (unsigned)conf.num_comp != conf.fact.size()) {
    std::cerr << "Invalid number of generator components in 'gen' tag.\n";
    return 1;
  }

  return 0;
}

//==============================================================================
//===== Dispatching reused tags to the right functions
//==============================================================================


template<typename Conf>
int tryTest(tinyxml2::XMLNode* current, Conf& conf) {
  //tinyxml2::XMLNode* node;
  if (!strcmp(current->Value(), "spectral")) {
    conf.criterion = LatticeTester::SPECTRAL;
    readSpectral(current, conf);
  } else if (!strcmp(current->Value(), "length")) {
    conf.criterion = LatticeTester::LENGTH;
    readLength(current, conf);
    std::cout << current->Value() << "\n";
  } else if (!strcmp(current->Value(), "beyer")) {
    std::cout << current->Value() << "\n";
  } else return 1;
  return 0;
}

//==============================================================================
//===== Different uses reding functions
//==============================================================================

template<typename Conf>
int readMK(tinyxml2::XMLNode* current, Conf& conf) {
  //if (current->NoChild()) return;
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!strcmp(node->Value(), "power")) {
      conf.power = true;
      conf.c1 = node->ToElement()->FirstAttribute()->IntValue();
    } else if (!strcmp(node->Value(), "range")) {
      conf.power = false;
      auto attr = node->ToElement()->FirstAttribute();
      conf.c1 = attr->IntValue();
      attr = attr->Next();
      conf.c2 = attr->IntValue();
    } else if (!strcmp(node->Value(), "k")) {
      conf.k = node->ToElement()->FirstAttribute()->IntValue();
    } else if (!strcmp(node->Value(), "safe")) {
      conf.safe = node->ToElement()->FirstAttribute()->BoolValue();
    } else if (!strcmp(node->Value(), "factor")) {
      conf.facto = node->ToElement()->FirstAttribute()->BoolValue();
    } else {
      std::cerr << "Invalid mk configuration.\n";
      return 1;
    }
    node = node->NextSibling();
  }
  return 0;
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
  tinyxml2::XMLNode* node = current->FirstChild();
  while (node) {
    if (!conf.gen_set && !strcmp(node->Value(), "gen")) {
      if (!readGenerator(node, conf)) conf.gen_set = true;
    } else if (!conf.test_set && !tryTest(node, conf)) {
      conf.test_set = true;
    } else if (!conf.proj_set && !strcmp(node->Value(), "proj")) {
      if (!readProj(node, conf)) conf.proj_set = true;
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
    if (!conf.gen_set && !strcmp(node->Value(), "gen")) {
      if (!readGenerator(node, conf)) conf.gen_set = true;
    } else if (!conf.test_set && !tryTest(node, conf)) {
      conf.test_set = true;
    } else if (!conf.proj_set && !strcmp(node->Value(), "proj")) {
      if (!readProj(node, conf)) conf.proj_set = true;
    } else if (!strcmp(node->Value(), "time")) {
      conf.timeLimit = node->ToElement()->FirstAttribute()->DoubleValue();
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
    if (types == "ZD" || types == "ZR") {
      MKSearch<NTL::ZZ> prog;
      if (readMK(current, prog)) return 1;
      return prog.FindMK();
    } else if (types == "LD") {
      MKSearch<std::int64_t> prog;
      readMK(current, prog);
      return prog.FindMK();
    }
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
    if (types == "LD") {
      LatTest<std::int64_t, double> prog;
      readTest(current, prog.conf);
      return prog.TestLat();
    } else if (types == "ZD") {
      LatTest<NTL::ZZ, double> prog;
      readTest(current, prog.conf);
      return prog.TestLat();
    } else if (types == "ZR") {
      LatTest<NTL::ZZ, NTL::RR> prog;
      readTest(current, prog.conf);
      return prog.TestLat();
    }
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

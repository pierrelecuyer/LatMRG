#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

#include "tinyxml2.h"

//#include "TestLattice.h"

#define MRGLATTICE_MAIN_EXEC
#include "ExecCommon.h"

#include "Seek.h"
#include "LatTest.h"
#include "FindMK.h"
#include "Period.h"

using namespace LatMRG;

// ===== Whole program configuration options ===================================

// Output stuff
std::ostream* out(&std::cout);
std::ofstream fout;
// Printing time or not. This is used to compare outputs
bool print_time = true;

std::string mode;

// Prints the program usage
void print_help() {
  *out << "Usage: MRGLattice [types] file1 file2 ...\n";
  *out << "types:\n"
    "  LD: use long and double types.\n"
    "  ZD: use NTL::ZZ and double types.\n"
    "  ZR: use NTL::ZZ and NTL::RR types.\n"
    "  This field can be ommited.\n";
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

/**
 * Fills the matrix `mat` with the values of the attributes of current.
 * This basically reshapes `mat` to a `dim*dim` size matrix and calls
 * `toVectString()` on the successive attributes of current. This will fail if
 * `current` has less than `dim` attributes or if any of those attributes does
 * not contain a vector of `dim` components.
 * */
template<typename Mat>
int toMatTag(tinyxml2::XMLElement* current, Mat& mat, int dim) {
  mat.SetDims(dim, dim);
  auto node = current->FirstAttribute();
  for (int i = 0; i < dim; i++) {
    if (node) {
      if (toVectString(node->Value(), mat[i], dim)) {
        std::cerr << "Cannot initialize line " << i+1 << " of the matrix.\n";
        return 1;
      }
    } else {
      std::cerr << "'matrix' tag expects " << dim << " attributes, but cannot "
        "find attribute " << i+1 << ".\n";
      return 1;
    }
    node = node->Next();
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

  int exponent, order;
  typename Conf::Int /*modulo,*/ basis, rest;

  std::string filem1 = "./tempm1" + i + std::to_string(rand());
  std::string filer = "./tempr" + i + std::to_string(rand());
  DecompType decompm1 = NO_DECOMP, decompr = NO_DECOMP;

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

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(order, node->FirstAttribute()->IntValue());
  } else {
    NTL::conv(order, 1);
  }

  if (mode == "seek") {
    node = current->FirstChildElement("method");
    if (node) {
      auto value = node->FirstAttribute();
      if (value) {
        if (!strcmp(value->Name(), "random")) {
          conf.search_mode[i] = "random";
          conf.coeff[i].SetLength(order);
          if (toVectString(value->Value(), conf.coeff[i], order)) return 1;
        } else if (!strcmp(value->Name(), "pow2")) {
          conf.search_mode[i] = "pow2";
          conf.coeff[i].SetLength(order+1);
          std::vector<long> tmp;
          tmp.resize(order+1);
          if (toVectString(value->Value(), tmp, order+1)) return 1;
          for (int j = 0; j < order; j++) conf.coeff[i][j] = typename Conf::Int(tmp[j+1]);
          conf.coeff[i][order] = typename Conf::Int(tmp[0]);
        } else if (!strcmp(value->Name(), "exhaust")) {
          conf.search_mode[i] = "exhaust";
          conf.coeff[i].SetLength(order);
          if (toVectString(value->Value(), conf.coeff[i], order)) return 1;
        } else {
          std::cerr << "Invalid attribute in tag 'method'.\n";
          return 1;
        }
      } else {
        std::cerr << "Tag 'method' has no attribute.\n";
        return 1;
      }
    } else {
      std::cerr << "No way to set coefficients in 'mrg' or 'lcg' tag. Add 'method'"
       " tag.\n";
      return 1;
    }
  }
  if (mode == "lattest") {
    node = current->FirstChildElement("coeff");
    if (node) {
      conf.coeff[i].SetLength(order);
      toVectString(node->FirstAttribute()->Value(), conf.coeff[i], order);
    } else {
      std::cerr << "No way to set coefficients in 'mrg' or 'lcg' tag. Add 'method' or"
        " 'coeff' tag.\n";
      return 1;
    }
  }

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf, i, filem1, filer, decompm1, decompr))
    std::cerr << "Non critical error in 'period' tag of 'mrg' tag.\n";

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
  NTL::matrix<typename Conf::Int> temp_mat;

  std::string filem1 = "./tempm1" + i + std::to_string(rand());
  std::string filer = "./tempr" + i + std::to_string(rand());
  DecompType decompm1 = NO_DECOMP, decompr = NO_DECOMP;

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

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(order, node->FirstAttribute()->Value());
  } else {
    std::cerr << "No 'order' tag in 'mrg' tag.\n";
    return 1;
  }

  if (mode == "seek") {
    std::cerr << "Currently no way to build matrix for MMRG implemented.\n";
    return 1;
  }
  if (mode == "lattest") {
    node = current->FirstChildElement("matrix");
    if (node) {
      temp_mat.SetDims(order, order);
      toMatTag(node, temp_mat, order);
    } else {
      std::cerr << "No way to set matrix in 'mmrg' tag. Add 'matrix' tag.\n";
      return 1;
    }
  }

  node = current->FirstChildElement("period");
  if (node && readGenPer(node, conf, i, filem1, filer, decompm1, decompr))
    std::cerr << "Non critical error in 'period' tag of 'mmrg' tag.\n";

  auto comp = new MRGComponent<typename Conf::Int>(basis, exponent, rest, order,
      decompm1, filem1.c_str(), decompr, filer.c_str());
  comp->setA(temp_mat);
  comp->set_type(MMRG);

  conf.fact.push_back(comp);

  return 0;
}

//==============================================================================

template<typename Conf>
int readMWC(tinyxml2::XMLNode* current, Conf& conf, int i) {
  tinyxml2::XMLElement* node;

  int exponent, order;
  typename Conf::Int modulo, basis, rest, m(0);

  std::string filem1 = "./tempm1" + i + std::to_string(rand());
  std::string filer = "./tempr" + i + std::to_string(rand());
  DecompType decompm1 = NO_DECOMP, decompr = NO_DECOMP;

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
    modulo = NTL::power(basis, exponent) + rest;
  } else {
    std::cerr << "No 'modulo' tag in 'mwc' tag.\n";
    return 1;
  }

  node = current->FirstChildElement("order");
  if (node) {
    NTL::conv(order, node->FirstAttribute()->IntValue());
  } else {
    std::cerr << "No 'order' tag in 'mwc' tag.\n";
    return 1;
  }

  if (mode == "seek") {
    node = current->FirstChildElement("method");
    if (node) {
      auto value = node->FirstAttribute();
      if (value) {
        if (!strcmp(value->Name(), "random")) {
          conf.coeff[i].SetLength(order);
          if (toVectString(value->Value(), conf.coeff[i], order)) return 1;
        }
      }
    } else {
      std::cerr << "Currently no way to search for coefficients for MWC.\n";
      return 1;
    }
  }
  if (mode == "lattest") {
    node = current->FirstChildElement("coeff");
    if (node) {
      conf.coeff[i].SetLength(order+1);
      toVectString(node->FirstAttribute()->Value(), conf.coeff[i], order+1);
    } else {
      std::cerr << "No way to set coefficients in 'mwc' tag. Add 'method' or"
        " 'coeff' tag.\n";
      return 1;
    }
  }

  if (mode == "lattest") {
    node = current->FirstChildElement("period");
    if (node && readGenPer(node, conf, i, filem1, filer, decompm1, decompr))
      std::cerr << "Non critical error in 'period' tag of 'mwc' tag.\n";

    for(int j = 0; j <= order; j++) {
      m += conf.coeff[i][j] * NTL::power(modulo, j);
    }
  }
  if (mode == "seek") {
    m = 2;
  }
  auto comp = new MRGComponent<typename Conf::Int>(m, order,
      decompm1, filem1.c_str(), decompr, filer.c_str());
  comp->set_type(MWC);
  comp->m_MWCb = modulo;
  comp->getB() = basis;
  comp->getE() = exponent;
  comp->getR() = rest;

  conf.fact.push_back(comp);

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
    if (LatticeTester::toPreRedString(conf.reduction, node->FirstAttribute()->Value())) std::cerr << "Invalid"
      "/inexistant attribute in 'reduction' tag.\n";
  }

  node = current->FirstChildElement("norma");
  if (node) {
    if (LatticeTester::toNormaString(conf.normaType, node->FirstAttribute()->Value())) std::cerr << "Invalid"
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
    if (LatticeTester::toPreRedString(conf.reduction, node->FirstAttribute()->Value())) std::cerr << "Invalid"
      "/inexistant attribute in 'reduction' tag.\n";
  }

  node = current->FirstChildElement("dual");
  if (node) {
    conf.use_dual = node->FirstAttribute()->BoolValue();
  }

  node = current->FirstChildElement("norm");
  if (node && LatticeTester::toNormString(conf.norm, node->FirstAttribute()->Value())) {
      std::cerr << "Ivalide first attribute value in node norm.\n";
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
    for (unsigned int i = 0; i < projDim.size(); i++) projDim[i] -= 1;
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
    decompm1 = DECOMP;
    decompr = DECOMP;
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
  auto tmp_node = current->FirstChildElement("numcomp");
  if (tmp_node) {
    conf.num_comp = tmp_node->FirstAttribute()->IntValue();
  }

  conf.coeff.resize(conf.num_comp);
  conf.period.resize(conf.num_comp);
  conf.search_mode.resize(conf.num_comp);

  srand(time(NULL));

  auto node = current->FirstChild();
  int count = 0;
  while (node && (count < conf.num_comp)) {
    if (!strcmp(node->Value(), "mrg")) {
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
      if (readMWC(node, conf, count)) {
        std::cerr << "'mwc' tag is not properly parametrized.\n";
        return 1;
      }
      count++;
    } else if (!strcmp(node->Value(), "lcg")) {
      if (readMRG(node, conf, count)) {
        std::cerr << "'lcg' tag is not properly parametrized.\n";
        return 1;
      }
      conf.fact[count]->set_type(LCG);
      count++;
    }

    node = node->NextSibling();
  }

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
      auto attr = node->ToElement()->FirstAttribute();
      conf.e = attr->IntValue();
      attr = attr->Next();
      conf.c1 = attr->IntValue();
    } else if (!strcmp(node->Value(), "range")) {
      conf.power = false;
      auto attr = node->ToElement()->FirstAttribute();
      conf.e = attr->IntValue();
      attr = attr->Next();
      conf.c1 = attr->IntValue();
      attr = attr->Next();
      conf.c2 = attr->IntValue();
    } else if (!strcmp(node->Value(), "k")) {
      conf.k = node->ToElement()->FirstAttribute()->IntValue();
    } else if (!strcmp(node->Value(), "safe")) {
      conf.safe = node->ToElement()->FirstAttribute()->BoolValue();
    } else if (!strcmp(node->Value(), "factor")) {
      conf.facto = node->ToElement()->FirstAttribute()->BoolValue();
    }
    node = node->NextSibling();
  }
  return 0;
}

//==============================================================================

template<typename Conf>
int readPeriod(tinyxml2::XMLNode* current, Conf& conf) {
  srand(time(NULL));
  conf.filem1 = "./tempm1" + std::to_string(rand());
  conf.filer = "./tempr" + std::to_string(rand());
  conf.decompm1 = DECOMP;
  conf.decompr = DECOMP;

  // Reading Modulo
  auto node = current->FirstChildElement("modulo");
  if (node) {
    if (!node->Attribute("basis")) {
      std::cerr << "No `basis` attribute in 'modulo' tag.\n";
      return 1;
    }
    NTL::conv(conf.m_b, node->Attribute("basis"));
    if (node->Attribute("exponent")) NTL::conv(conf.m_e, node->Attribute("exponent"));
    else conf.m_e = 1;
    if (node->Attribute("rest")) NTL::conv(conf.m_r, node->Attribute("rest"));
    else conf.m_r = 0;
    conf.m_m = NTL::power(conf.m_b, conf.m_e) + conf.m_r;
  } else {
    std::cerr << "No 'modulo' tag in 'period' tag.\n";
    return 1;
  }

  // Reading type
  node = current->FirstChildElement("type");
  if (node && toGenString(conf.type, node->FirstAttribute()->Value())) {
    std::cerr << "Invalid attribute " << node->FirstAttribute()->Name()
      << " value in 'type' tag.\n";
  }

  // Reading order
  if (conf.type != LCG) {
    node = current->FirstChildElement("order");
    if (node) {
      NTL::conv(conf.m_k, node->FirstAttribute()->IntValue());
      conf.m_a.SetLength(conf.m_k+1);
    } else {
      std::cerr << "No 'order' tag in 'period' tag.\n";
      return 1;
    }
  } else {
    conf.m_k = 1;
  }

  if (conf.type == MRG) {
    NTL::vector<typename Conf::Int> a;
    a.SetLength(conf.m_k);
    node = current->FirstChildElement("mult");
    if (!node) {
      std::cerr << "No 'mult' tag. Cannot check for full period without"
        " coefficients.\n";
      return 1;
    }
    if (toVectString(node->FirstAttribute()->Value(), a, conf.m_k)) {
      std::cerr << "Error in 'mult' tag.\n";
    } else {
      for (int i = 0; i < conf.m_k; i++) conf.m_a[i+1] = a[i];
    }
  } else if (conf.type == MMRG) {
    node = current->FirstChildElement("matrix");
    if (!node) {
      std::cerr << "No 'matrix' tag. Cannot check for full period without"
        " multiplier matrix.\n";
      return 1;
    }
    if (toMatTag(node, conf.m_A, conf.m_k)) {
      std::cerr << "Error reading 'matrix' tag.\n";
      return 1;
    }
  } else if (conf.type == MWC) {
    node = current->FirstChildElement("mult");
    if (!node) {
      std::cerr << "No 'mult' tag. Cannot check for full period without"
        " coefficients.\n";
      return 1;
    }
    conf.m_a.SetLength(conf.m_k+1);
    if (toVectString(node->FirstAttribute()->Value(), conf.m_a, conf.m_k+1))
      std::cerr << "Error in 'mult' tag.\n";
  } else if (conf.type == LCG) {
    node = current->FirstChildElement("mult");
    if (!node) {
      std::cerr << "No 'mult' tag. Cannot check for full period without"
        " coefficient.\n";
      return 1;
    }
    NTL::conv(conf.m_mult, node->FirstAttribute()->Value());
  }

  // Read decomposition files configuration
  int err = 0;
  node = current->FirstChildElement("m1");
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
  if (err) return 1;

  return 0;
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
    } else if (!strcmp(node->Value(), "detail")) {
      conf.detail = node->ToElement()->FirstAttribute()->IntValue();
    }
    node = node->NextSibling();
  }
}


//==============================================================================
//===== Main reading function
//==============================================================================

template<typename Int, typename Dbl>
int readFile(const char* filename) {
  mode = "";

  tinyxml2::XMLDocument doc;
  doc.LoadFile(filename);
  tinyxml2::XMLNode* current;
  current = doc.FirstChild();
  out = &std::cout;
  if (doc.FirstChildElement("out")) {
    auto attr = doc.FirstChildElement("out")->FirstAttribute();
    if (attr) {
      fout = std::ofstream(attr->Value());
      out = &fout;
    } else {
      std::string name = filename;
      while (name.back() != '.') name.pop_back();
      name += "res";
      fout = std::ofstream(name);
      out = &fout;
    }
  }
  print_time = true;
  if (doc.FirstChildElement("print_time")) {
    auto attr = doc.FirstChildElement("print_time")->FirstAttribute();
    if (attr) print_time = attr->BoolValue();
  }
  while (current) {
    if (!strcmp(current->Value(), "mk")) {
      mode = "mk";
      MKSearch<Int> prog;
      if (readMK(current, prog)) return 1;
      return prog.FindMK();
    } else if (!strcmp(current->Value(), "period")) {
      mode = "period";
      MaxPeriod<Int, Dbl> prog;
      if (readPeriod(current, prog)) return 1;
      return prog.CheckPeriod();
    } else if (!strcmp(current->Value(), "seek")) {
      mode = "seek";
      ConfigSeek<Int, Dbl> conf;
      readSeek(current, conf);
      if (conf.num_comp > 1) {
        SeekMain<ComboLattice<Int, Dbl>> prog(conf);
        return prog.Seek(LatMRGSeek::ComboLatticeFinder<Int, Dbl>::getFunction);
      } else if (conf.fact[0]->get_type() == MRG) {
        SeekMain<MRGLattice<Int, Dbl>> prog(conf);
        if (search_mode == "exhaust") {
          return prog.Seek(LatMRGSeek::MRGLatticeExhaust<Int, Dbl>::nextGenerator);
        } else if (search_mode == "pow2") {
          return prog.Seek(LatMRGSeek::nextGeneratorPow2);
        }
        return prog.Seek(LatMRGSeek::nextGenerator);
      } else if (conf.fact[0]->get_type() == LCG) {
        SeekMain<MRGLattice<Int, Dbl>> prog(conf);
        return prog.Seek(LatMRGSeek::nextGenerator);
      } else if (conf.fact[0]->get_type() == MWC) {
        SeekMain<MWCLattice<Int, Dbl>> prog(conf);
        return prog.Seek(LatMRGSeek::nextGenerator);
      } else if (conf.fact[0]->get_type() == MMRG) {
        SeekMain<MMRGLattice<Int, Dbl>> prog(conf);
        return prog.Seek(LatMRGSeek::nextGenerator);
      }
      std::cerr << "Seek exited prematurely.\n";
      return 1;
    } else if (!strcmp(current->Value(), "lattest")) {
      mode = "lattest";
      LatTest<Int, Dbl> prog;
      readTest(current, prog.conf);
      return prog.TestLat();
    }
    current = current->NextSibling();
  }
  return 0;
}

//==============================================================================

int main(int argc, char** argv) {
  if (argc < 2) {
    print_help();
    return 1;
  }

  // Check types
  if (!(strcmp(argv[1],"ZD") && strcmp(argv[1], "LD") && strcmp(argv[1], "ZR"))) types = argv[1];
  else if (readFile<NTL::ZZ, double>(argv[1])) {
    std::cerr << "File " << argv[1] << " exited with errors.\n";
  }

  // Run all specified files
  for (int i = 2; i < argc; i++) {
    if (types == "LD") {
      if (readFile<std::int64_t, double>(argv[i])) {
        std::cerr << "File " << argv[i] << " exited with errors.\n";
      }
    } else if (types == "ZD") {
      if (readFile<NTL::ZZ, double>(argv[i])) {
        std::cerr << "File " << argv[i] << " exited with errors.\n";
      }
    } else if (types == "ZR") {
      if (readFile<NTL::ZZ, NTL::RR>(argv[i])) {
        std::cerr << "File " << argv[i] << " exited with errors.\n";
      }
    }
  }
  return 0;
}

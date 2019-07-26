/**
 * Parameters list:
 * Lat
 * GenType     (MRG)
 * Criterion   (SPECTRAL)
 * UseDual     (true)
 * Reduction   (FULL)
 * NormaType   (NONE)
 * Projections (no default)
 * timeLimit   (3600s)
 * detail      (medium)
 * period      (no)
 * generator detail (no default)
 * MRG
 * a
 * m
 * k
 * MWC
 * a
 * m
 * k
 * COMBO
 * J
 * a_j
 * m_j
 * k_j
 * MMRG
 * A
 * m
 * k
 *
 * Seek
 *
 * Period
 *
 * FindMK
 * */
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>

#include "tinyxml2.h"

#include "TestLattice.h"

#define MRGLATTICE_MAIN_EXEC
#include "ExecCommon.h"

//using namespace LatMRG;

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

int parse_option(char* option, int& i, int argc) {
  char name[40];
  char argument[40] = "";
  // Splitting the option name and the argument
  if (option[1] == '-') {
    strcpy(name, option+2);
    for (unsigned int i = 0; i < strlen(name); i++) {
      if (name[i] == '=') {
        strcpy(argument, name+i+1);
        name[i] = '\0';
        break;
      }
    }
  } else {
    strcpy(name, option+1);
    if (name[1] == '=') strcpy(argument, name+2);
    name[1] = '\0';
  }

  if (!(strcmp(name, "h") && strcmp(name, "help"))) {
    return 1;
  } else if (!(strcmp(name, "o") && strcmp(name, "output"))){
    if (strlen(argument) == 0) {
      std::cerr << "Missing argument to option -o, --output\n";
      return 1;
    }
    fout.open(argument);
    out = &fout;
    return 0;
  } else if (!(strcmp(name, "v") && strcmp(name, "verbose"))) {
    std::cerr << "No implementation for verbose option.\n";
    return 0;
  } else if (!(strcmp(name, "t") && strcmp(name, "types"))) {
    if (strlen(argument) == 0) {
      std::cerr << "Missing argument to option -t, --types\n";
      return 1;
    }
    std::cerr << "No implementation for types option.\n";
    return 0;
  } else {
    std::cerr << "Unknown option\n";
  }
  return 1;
}

int parse_param(char* option, int& i, int argc) {
  char name[40];
  char argument[40] = "";
  // Splitting the option name and the argument
  if (option[1] == '-') {
    strcpy(name, option+2);
    for (unsigned int i = 0; i < strlen(name); i++) {
      if (name[i] == '=') {
        strcpy(argument, name+i+1);
        name[i] = '\0';
        break;
      }
    }
  } else {
    strcpy(name, option+1);
    if (name[1] == '=') strcpy(argument, name+2);
    name[1] = '\0';
  }

  // The first layer of if/elses filters between the different modes
  if (exec_mode == "lat") {
    if (!strcmp(name, "gentype")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --gentype\n";
        return 1;
      }
      return toGenString(gen_type, std::string(argument));
    } else if (!strcmp(name, "criterion")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --criterion\n";
        return 1;
      }
      return toCriterionString(crit_type, std::string(argument))
    } else if (!strcmp(name, "dual")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --dual\n";
        return 1;
      }
      if (!strcmp(argument, "true")) {
        dual = true;
        return 0;
      } else if (!strcmp(argument, "false")){
        dual = false;
        return 0;
      }
    } else if (!strcmp(name, "reduction")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --reduction\n";
        return 1;
      }
      return toRedString(red_type, std::string(argument));
    } else if (!strcmp(name, "normalizer")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --normalizer\n";
        return 1;
      }
    } else if (!strcmp(name, "time")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --time\n";
        return 1;
      }
    } else if (!strcmp(name, "projections")) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter --projections\n";
        return 1;
      }
    } else if (!(strcmp(name, "a") && strcmp(name, "vector"))) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter -a --vector\n";
        return 1;
      }
    } else if (!(strcmp(name, "m") && strcmp(name, "modulo"))) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter -m --modulo\n";
        return 1;
      }
    } else if (!(strcmp(name, "k") && strcmp(name, "order"))) {
      if (strlen(argument) == 0) {
        std::cerr << "Missing argument to parameter -k --order\n";
        return 1;
      }
    }
  } else if (exec_mode == "seek") {
  } else if (exec_mode == "mk") {
  } else if (exec_mode == "period"){
  }
  return 1;
}

int parse_mode(char* mode) {
  if (strcmp(mode, "lat")) {
    exec_mode = "lat";
    options = false;
    return 0;
  } else if (strcmp(mode, "seek")) {
    exec_mode = "seek";
    options = false;
    return 0;
  } else if (strcmp(mode, "mk")) {
    exec_mode = "mk";
    options = false;
    return 0;
  } else if (strcmp(mode, "period")) {
    exec_mode = "period";
    options = false;
    return 0;
  } else {
    std::cerr << "Invalid mode\n";
  }
  return 1;
}

int parse_arguments(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Missing arguments\n";
    return 1;
  }
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (options) {
        if (parse_option(argv[i], i, argc)) return 1;
      } else {
        if(parse_param(argv[i], i, argc)) return 1;
      }
    } else {
      if(parse_mode(argv[i])) return 1;
    }
  }
  return 0;
}

int main(int argc, char** argv) {
  if (parse_arguments(argc, argv)) {
    print_help();
  }
  if (exec_mode == "lat") {
    return TestLattice(gen_type, crit_type, dual, red_type, norma_type,
        time_limit, detail, period, *proj, gen_string);
  } else if (exec_mode == "seek") {
  } else if (exec_mode == "mk") {
  } else if (exec_mode == "period") {
  }
  return 1;
}

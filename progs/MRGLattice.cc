#include <cstring>
#include <string>
#include <iostream>
#include <fstream>

//using namespace LatMRG;

bool options = true;
std::string exec_mode;

// Output stuff
std::ofstream fout;
std::ostream* out(&std::cout);

// Prints the program usage
void print_help() {
  std::cout << "Usage: MRGLattice [OPTIONS] MODE PARAMETERS\n\n";
  std::cout << "Available options:\n\t" << "-h, --help\t\t" << "Print this"
   " message.\n" << "\t-o, --output=FILE \t" << "Output to FILE\n" << "\n";
  std::cout << "Available modes:\n\t" << "lat\t\t\t" << "Study the"
    " lattice of a MRG generator.\n\t" << "seek\t\t\t" << "Search for"
    " new MRG generators.\n\t" << "mk\t\t\t" << "Find compatible "
    "modulus and order for MRG generators.\n\t" << "period\t\t\t" <<
    "Test the period of a MRG generator.\n";
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
      std::cout << "Missing argument to option -o, --output\n";
      return 1;
    }
    fout.open(argument);
    out = &fout;
    return 0;
  } else {
    std::cout << "Unknown option\n";
  }
  return 1;
}

int parse_param(char* option, int& i, int argc) {
  return 0;
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
    std::cout << "Invalid mode\n";
  }
  return 1;
}

int parse_arguments(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Missing arguments\n";
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
  *out << "hi\n";
  return 0;
}

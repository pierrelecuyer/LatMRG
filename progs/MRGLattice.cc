#include <iostream>

//using namespace LatMRG;

// Prints the program usage
void print_help() {
  std::cout << "Usage: MRGLattice MODE PARAMETERS [OPTIONS]\n\n";
  std::cout << "Available modes:\n\t" << "lat\t\t" << "Study the"
    " lattice of a MRG generator.\n\t" << "seek\t\t" << "Search for"
    " new MRG generators.\n\t" << "mk\t\t" << "Find compatible "
    "modulus and order for MRG generators.\n\t" << "period\t\t" <<
    "Test the period of a MRG generator.\n\t" << "help\t\t" <<
    "Print this message.\n";
}

int main(int argc, char** argv) {
  if (argc < 3) {
    print_help();
    return 0;
  }
  return 0;
}

/*
 * This program look for integers smaller to specific powers of two. Look at
 * the manual for more details.
 * */

#include "latmrg/Primes.h"
#include "latmrg/ParamReaderExt.h"

#include <NTL/ZZ.h>

#include <stdint.h>
#include <fstream>

using namespace std;
using namespace LatMRG;
using namespace NTL;

int main(int argc, char** argv)
{
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <I,Z> <data file>" << endl;
    cout << 1 <<endl;
    return 1;
  }
  string types(argv[1]);
  if (types != "I" && types != "Z") {
    cout << "Usage: " << argv[0] << " <I,Z> <data file>" << endl;
    cout << 2 <<endl;
    return 1;
  }
  string fname(argv[2]);
  ofstream fout (fname + ".res");
  fname += ".dat";

  ParamReaderExt<int64_t, double> reader(fname);
  reader.getLines();
  int k, e, c1, c2;
  bool power, safe;
  int ln = 0;
  reader.readBool(power, ln, 0);
  reader.readInt(k, ++ln, 0);
  reader.readInt(e, ++ln, 0);
  if (power) {
    reader.readInt(c1, ++ln, 0);
  } else {
    reader.readInt(c1, ++ln, 0);
    reader.readInt(c2, ln, 1);
  }
  reader.readBool(safe, ++ln, 0);

  //   primes.find (3, 39, 3, true, false, fout);

  // trouver 3 modules m proches de 2^31 pour des MRGs d'ordre k
  if (types == "I") {
    Primes<std::int64_t> primes;
    if (power) primes.find (k, e, c1, safe, false, fout);
    else primes.find (k, e, c1, c2, safe, false, fout);
  } else if (types == "Z") {
    Primes<NTL::ZZ> primes;
    if (power) primes.find (k, e, c1, safe, false, fout);
    else primes.find (k, e, c1, c2, safe, false, fout);
  }
  return 0;
}

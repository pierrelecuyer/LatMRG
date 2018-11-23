#include "latmrg/Primes.h"
#include "latmrg/ParamReaderExt.h"

#include "latticetester/Types.h"

#include <string>

using namespace std;
using namespace LatMRG;

//===========================================================================

namespace
{

  typedef struct {
    int k;
    int e;
    long c1;
    long c2;
    bool safe;
    bool fac;
  } ParamType;


  //===========================================================================

  template<typename Int, typename Dbl>
    void read (ParamReaderExt<Int, Dbl> & reader, ParamType & par)
    {
      reader.getLines ();
      int ln = 0;
      reader.readInt (par.k, ++ln, 1);
      reader.readInt (par.e, ++ln, 1);
      reader.readLong (par.c1, ++ln, 1);
      reader.readLong (par.c2, ++ln, 1);
      reader.readBool (par.safe, ++ln, 1);
      reader.readBool (par.fac, ++ln, 1);
    }

  void ID(string fname) {
    fname += ".dat";
    ParamReaderExt<std::int64_t, double> reader (fname);
    ParamType par;
    read (reader, par);

    string oname(fname);
    oname += ".res";

    ofstream fout (oname.c_str());
    Primes<std::int64_t> primes;
    primes.find (par.k, par.e, par.c1, par.c2, par.safe, par.fac, fout);
  }

  void ZD(string fname) {
    fname += ".dat";
    ParamReaderExt<NTL::ZZ, double> reader (fname);
    ParamType par;
    read (reader, par);

    string oname(fname);
    oname += ".res";

    ofstream fout (oname.c_str());
    Primes<NTL::ZZ> primes;
    primes.find (par.k, par.e, par.c1, par.c2, par.safe, par.fac, fout);
  }

  void ZR(string fname) {
    fname += ".dat";
    ParamReaderExt<NTL::ZZ, NTL::RR> reader (fname);
    ParamType par;
    read (reader, par);

    string oname(fname);
    oname += ".res";

    ofstream fout (oname.c_str());
    Primes<NTL::ZZ> primes;
    primes.find (par.k, par.e, par.c1, par.c2, par.safe, par.fac, fout);
  }

}  // namespace


//===========================================================================

int main(int argc, char** argv)
{
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <data file>" << endl;
    return 1;
  }
  string types(argv[1]);
  string fname(argv[2]);
  if (types == "ID") {
    ID(fname);
  } else if (types == "ZD") {
    ZD(fname);
  } else if (types == "ZR") {
    ZR(fname);
  }
  return 0;
}

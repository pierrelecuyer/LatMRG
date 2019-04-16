#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include "latticetester/Types.h"
#include "latmrg/LatTestAll.h"


using namespace std;
using namespace LatMRG;


//==========================================================================

int main (int argc, char **argv)
{

  struct stat buf;    // properties of a file or directory
  LatTestAll<std::int64_t, double> testallID;
  LatTestAll<NTL::ZZ, double> testallZD;
  LatTestAll<NTL::ZZ, NTL::RR> testallZR;
  int status = 0;
  string types(argv[1]);

  if (argc < 3 || (types != "ID" && types != "ZD" && types != "ZR")) {
    cerr << "\n*** Usage:\n   "
      << argv[0] << " <ID,ZD,ZR> data_file1 data_file2 ...." << endl
      << "or\n   "
      << argv[0] << " <ID,ZD,ZR> dir1 dir2 ...." << endl
      << endl;
    return -1;
  }

  for (int j = 2; j < argc; j++) {
    // Do the test for each data file or directory on the command line
    stat(argv[j], &buf);
    if (0 != S_ISDIR(buf.st_mode)) {
      std::cout << "Is dir\n";
      if (types == "ID") {
        status |= testallID.doTestDir (argv[j]);
      } else if (types == "ZD") {
        status |= testallZD.doTestDir (argv[j]);
      } else if (types == "ZR") {
        status |= testallZR.doTestDir (argv[j]);
      }
    } else {
      string dataname(argv[j]);
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      if (0 != S_ISREG(buf.st_mode)){    // data file
        if (types == "ID") {
          status |= testallID.doTest (argv[j]);
          std::cout << "end\n";
        } else if (types == "ZD") {
          int temp = testallZD.doTest(argv[j]);
          std::cout << "end\n";
          status |= temp;
        } else if (types == "ZR") {
          status |= testallZR.doTest (argv[j]);
          std::cout << "end\n";
        }
      }
    }
  }

  return status;
}

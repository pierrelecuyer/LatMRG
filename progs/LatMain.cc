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

  if (argc < 3) {
    cerr << "\n*** Usage:\n   "
      << argv[0] << " types data_file1 data_file2 ...." << endl
      << "or\n   "
      << argv[0] << " types dir1 dir2 ...." << endl
      << endl;
    return -1;
  }

  struct stat buf;    // properties of a file or directory
  LatTestAll<std::int64_t, double> testallID;
  LatTestAll<NTL::ZZ, double> testallZD;
  LatTestAll<NTL::ZZ, NTL::RR> testallZR;
  int status = 0;
  string types(argv[1]);
  std::cout << types << std::endl;

  //char *testfile = "/Users/Erwan1/projects/github/LatMRG/latZZDD_test2";
  //status |= testall.doTest (testfile);

  for (int j = 2; j < argc; j++) {
    // Do the test for each data file or directory on the command line
    stat(argv[j], &buf);
    if (0 != S_ISDIR(buf.st_mode)) {
      if (types == "ID") {
        status |= testallID.doTestDir (argv[j]);
      } else if (types == "ZD") {
        status |= testallZD.doTestDir (argv[j]);
      } else if (types == "ZR") {
        status |= testallZR.doTestDir (argv[j]);
      }
      else {
        string dataname(argv[j]);
        dataname.append(".dat");
        stat(dataname.c_str(), &buf);
        if (0 != S_ISREG(buf.st_mode)){    // data file
          if (types == "ID") {
            status |= testallID.doTestDir (argv[j]);
          } else if (types == "ZD") {
            status |= testallZD.doTestDir (argv[j]);
          } else if (types == "ZR") {
            status |= testallZR.doTestDir (argv[j]);
          }
        }
      }
    }
  }

  return status;
}

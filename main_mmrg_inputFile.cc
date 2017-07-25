//
//  main_for_mmrg_input.cc
//  main program to test LatMRG execution
//  on an external input file for MMRG
//

// Include Header
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

// Include LatticeTester Header
/*
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"
#include "latticetester/Types.h"
#include "latticetester/ParamReader.h"
#include "latticetester/LatticeTesterConfig.h"
#include "latticetester/LatticeAnalysis.h"
*/

// Include LatMRG Header
#include "latmrg/LatTestAll.h"

// Include NTL Header
#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

// Include Boost Header
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

// Include Random Generator of MRG Matrix an tools
//#include "SimpleMRG.h"
//#include "Tools.h"

using namespace std;
using namespace NTL;
using namespace LatMRG;


//==========================================================================

int main ()
{
   // Erwan
   //string testLocation = "ton_path_vers_dossier_input_files_latmrg";
   //string testLocation = "ton_path_vers_dossier_input_files_latmrg/input_file_name";
   
   // Paul
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles";
   
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/latZZDD_test1";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/latZZDD_test2";
   
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/lacunaryMRG_test1";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/lacunaryMRG_test2";
   
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test1";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test2";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test2_bis";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test3";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test4";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test5";
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test6";
   string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mmrg_test7";
   
   //string testLocation = "/Users/paulwambergue/UdeM/latmrg/inputTestFiles/mrg_order2_test1";
   
   
   struct stat buf; // properties of a file or directory
   LatTestAll testall;
   int status = 0;
   
   stat(testLocation.c_str(), &buf);
   
   if (0 != S_ISDIR(buf.st_mode)) //directory
      status |= testall.doTestDir (testLocation.c_str());
   else { //file
      string dataname(testLocation.c_str());
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      status |= testall.doTest (testLocation.c_str());
   }
   
   return status;
}






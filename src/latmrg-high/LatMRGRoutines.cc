// This file is part of LatMRG.
//
// LatMRG
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <fnmatch.h>
#include <dirent.h>
#include <sys/types.h>
#include <cerrno>
#include <vector>
#include <string>
#include <iostream>

#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/WriterRes.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaPalpha.h"

#include "latmrg/LatConfig.h"
#include "latmrg/ParamReaderLat.h"
#include "latmrg/KorobovLattice.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/LatTestAll.h"
#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/Formatter.h"
#include "latmrg/ReportHeaderLat.h"
#include "latmrg/ReportFooterLat.h"
#include "latmrg/ReportLat.h"
#include "latmrg/TestProjections.h"
#include "latmrg-high/LatMRGRoutines.h"

using namespace std;
// using namespace NTL;
using namespace LatMRG;
using namespace LatticeTester;

namespace LatMRG {

  void printResult(const std::vector<double> & result, const int & fromDim)
  {
    cout << "Result: " << endl;
    for (unsigned int i = 0; i < result.size(); ++i) {
      cout << "  dim " << fromDim+i << " = " << sqrt(result[i]) << endl;
    }

  }
}

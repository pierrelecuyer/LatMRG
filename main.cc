//
//  main.cpp
//  LatMRG
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//


// select pre compiling options
//----------------------------------------------------------------------------------------

#define PRINT_CONSOLE

//----------------------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <num.h>

#include "latmrg/IntLattice.h"
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

#include "latmrg/LatTestAll.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>



#include "SimpleMRG.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;


int main (int argc, char *argv[])
{
   
   
   LatMRG::IntLattice *lattice = 0;
   
   BMat A;
   A.resize(2, 2);
   A(0,0) = 1;
   A(1,0) = 0;
   A(0,1) = 33;
   A(1,1) = 89;
   
   BMat B;
   B.resize(2, 2);
   B(0,0) = 0;
   B(1,0) = 0;
   B(0,1) = 0;
   B(1,1) = 0;
   
   MScal m(89);
   
   CalcDual(A, B, 2, m);
   
   cout << "A : " << endl;
   cout << A << endl;
   cout << "\n \n B : " << endl;
   cout << B << endl;
   
   
   
   cout << "Hello World" << endl;

   return 0;
}


#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#define LATMRG_SEEK

#include <iostream>
#include <cstdint>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/ConfigSeek.h"
#include "latmrg/Seek.h"

/**
 * This example implements a Seek for MRGs. It is far from being
 * the final version but should only be seen as a starting point
 */

using namespace LatMRG;



// Everything which is needed for the configuration is set first
// Only the last two lines are really performing the seek.
int main() {
  cout << "The purpose is to test Seek.h.\n";
  ConfigSeek<Int, Real> conf;
  NTL::Vec<int64_t> t; // The t-vector for the FOM.
  t.SetLength(3);
  t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
  t[1] = 5;    // Then pairs and triples up to coord. 5.
  t[2] = 5;
  NTL::ZZ m;
  m = 9223372036854773561;
  conf.configFOM.norma = new NormaBestLat(log(m), 1, 16);
  conf.configFOM.t = t;
  conf.genType = MRG;
  conf.numComp = 1;
  conf.genComponents.resize(conf.numComp);
  conf.genComponents[0] = new ConfigSeekMRG<Int, Real>();
  asMRG(conf.genComponents[0])->modulus = m;  
  NTL::Vec<NTL::ZZ> b;
  b.SetLength(3);
  b[1] = 1145902849652723;  
  b[2] = 5189485190151516; 
  NTL::Vec<NTL::ZZ> c;
  c.SetLength(3);
  c[1] = 1145902849652725;
  c[2] = 5189485190151518; 
  asMRG(conf.genComponents[0])->lowBoundaries = b;
  asMRG(conf.genComponents[0])->highBoundaries = c;
  asMRG(conf.genComponents[0])->order = b.length() - 1;
  conf.configFOM.max_gen = 4;
  conf.maxdim = 16;
  Seek<MRGLattice<Int, double>> seeker(conf);
  seeker.PerformSeek();
  return(0);
}

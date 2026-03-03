#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"

/**
 * This example makes speed comparisons for factorizing large numbers 
 * with the external tools 'yafu' and 'msieve'. It analyzes how the
 * performace of the two tools compares for different sizes of numbers
 * ranging from 10 to 70 digits.
 */

using namespace LatticeTester;

const int64_t numTestCases = 7;
const int64_t numMeth = 2;
// The following numbers are all primes and serve as the starting point for performing factorization. They consist of 10, 20, 30, 40, 50, 60 and 70 digits. 
const string startingNumbers[numTestCases] = { "3333555227", "10089886811898868001", "113333555555777777799999799999", "2030507011013017019023029031037041043047",
                                              "22953686867719691230002707821868552601124472329079", "111111111111111112222333333333334445556677777777777889999999", 
                                              "2030507011013017019023029031037041043047053059061067071073079083089097"};
const int64_t noFactorizations[numTestCases] = { 100, 100, 100, 100, 10, 10, 10}; // Numbers of factorizations to be performed for each of the number of digits
const string names[numMeth] = { "msieve              ", "yafu                "};
std::chrono::time_point<std::chrono::steady_clock> startTime, endTime;
std::chrono::milliseconds elapsedTime[numMeth][numTestCases];

void printResults(int64_t numSizes) {
   int64_t d;
   std::cout << "Timings for 'msieve' and 'yafu' to perform the factorization of 'n' numbers with 'd' digits (in milliseconds). \n \n";
   std::cout << "Digits:               ";
   for (d = 0; d < numTestCases; d++)
      std::cout << std::setw(8) << startingNumbers[d].length() << "  ";
   std::cout << "\n";
   std::cout << "No factorizations:    ";
   for (d = 0; d < numTestCases; d++)
      std::cout << std::setw(8) << noFactorizations[d] << "  ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      for (d = 0; d < numTestCases; d++)
         std::cout << std::setw(9) << elapsedTime[meth][d].count() << " ";
      std::cout << "\n";
   }
   std::cout << "\n";
}

int main() {
    NTL::ZZ b;
    b = 1;
    LatMRG::IntFactorization<Int> fact(b);
    for (int i = 0; i < numTestCases; i++)
    {
      b = conv<ZZ>(startingNumbers[i].c_str());
      // Run msieve
      startTime = std::chrono::steady_clock::now();
      for (int j = 0; j < noFactorizations[i]; j++)
      {
         fact.setNumber(b+j);
         fact.factorize(false);
      }
      endTime = std::chrono::steady_clock::now();      
      elapsedTime[0][i] =  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
      // Run yafu
      startTime = std::chrono::steady_clock::now();
      for (int j = 0; j < noFactorizations[i]; j++)
      {
         fact.setNumber(b+j);
         fact.factorize(true);
      }
      endTime = std::chrono::steady_clock::now();      
      elapsedTime[1][i] =  std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);   
    }
  printResults(numTestCases);
  return(0);
}

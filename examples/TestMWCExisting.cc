/**
 * TestMWCExisting: Testing previously-proposed MWC generators.
 */
#define TYPES_CODE  ZX     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <string>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/Random.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"

using namespace LatticeTester;
using namespace LatMRG;

//int64_t maxdim = 48;
//int64_t maxdim2 = 52;
LCGLattice<Int, Real> lcg(52);
ReducerBB<Int, Real> red(lcg);
WeightsUniform weights(1.0);
NormaBestLat normaDual(48);
FigureOfMeritDualM<Int, Real> fomdual(weights, normaDual, &red);

int64_t maxdim1 = 30;
LCGLattice<Int, Real> lcg1(30, L1NORM);
ReducerBB<Int, Real> red1(lcg1);
NormaMinkL1 normaDual1(30);
FigureOfMeritDualM<Int, Real> fomdual1(weights, normaDual1, &red1);

// Tests a MWC with b = 2^e, order k, coefficients aa, in up to maxdim dimensions.
template<typename Int, typename Real>
static void testMWCProposal(std::string name, int64_t k, int64_t e, IntVec aa, int64_t maxdim) {
   Int b = NTL::power(Int(2), e);
   std::cout << "_____________________________________________________________________\n";
   std::cout << name << "\n\n";
   std::cout << "MWC of order k = " << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   std::cout << "Coefficients aa = " << aa << "\n";
   Int m = computeLCGModulusMWC(b, aa);
   std::cout << "Modulo m = " << m << "\n";
   std::cout << "Log_2(m) = " << Lg(m) << "\n";
   if (mIsPrime<Int>(m) <= 1) std::cout << " m is prime.\n";
   else std::cout << " m is NOT prime.\n";
   if (IntFactor<Int>::isPrime((m - Int(1)) / Int(2)) <= 1) std::cout << " (m=1)/2 = "
         << (m - 1) / 2 << " is prime.\n";
   else std::cout << " (m=1)/2 = " << (m - 1) / 2 << " is NOT prime.\n";
   Int suma2 = aa[0] * aa[0];
   for (int64_t j = 1; j <= k; j++)
      suma2 += aa[j] * aa[j];
   std::cout << "a_0^2 + ... + a_k^2 = " << conv<double>(suma2) << "\n";
   std::cout << "bound k dim: (b^2 + 1)^{-1/2} = " << 1.0 / sqrt(conv<double>(b * b + 1)) << "\n";
   std::cout << "bound k+1 dim: (a_0^2 + ... + a_k^2)^{-1/2} = " << 1.0 / sqrt(conv<double>(suma2))
         << "\n\n";
   std::cout << "spectral test in dual lattice, L2 norm:\n";

   // Apply spectral test.
   lcg.setModulus(m);
   lcg.seta(b, maxdim);
   fomdual.getNormalizer()->computeBounds(-log(m), 1);
   fomdual.computeMeritSucc(lcg, k + 1, maxdim);  // Start in k+1 dim.
   std::cout << "\n";

   if (k + 2 < maxdim1) {
      std::cout << "spectral test in dual lattice, L1 norm:\n";
      lcg1.setModulus(m);
      lcg1.seta(b, maxdim1);
      fomdual1.getNormalizer()->computeBounds(-log(m), 1);
      fomdual1.computeMeritSucc(lcg1, k + 1, k + 3);  // Start in k+1 dim.
      std::cout << "\n";
   }
}

int main() {
   std::cout << "\n=============================================================\n";
   std::cout << "TestMWCExisting: ";
   std::cout << "We test some previously-proposed MWC generators of order k and modulo b.\n";
   fomdual.setVerbosity(4, 0.1);
   fomdual1.setVerbosity(4, 0.1);
   IntVec aa;
   int64_t maxdim;

   // ==============================================================================

   std::cout << "\n=============================================================\n";
   std::cout << "MWC generators proposed by Marsaglia (1994) .\n\n";
   maxdim = 16;

   // Example 1 from Marsaglia (1994)
   aa.SetLength(9);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>(1941);
   aa[2] = conv<ZZ>(1860);
   aa[3] = conv<ZZ>(1812);
   aa[4] = conv<ZZ>(1776);
   aa[5] = conv<ZZ>(1492);
   aa[6] = conv<ZZ>(1215);
   aa[7] = conv<ZZ>(1066);
   aa[8] = conv<ZZ>(12013);
   testMWCProposal<Int, Real>("Marsa 16a", 8, 16, aa, maxdim);

   // Example 2 from Marsaglia (1994)
   aa.SetLength(9);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>(1111);
   aa[2] = conv<ZZ>(2222);
   aa[3] = conv<ZZ>(3333);
   aa[4] = conv<ZZ>(4444);
   aa[5] = conv<ZZ>(5555);
   aa[6] = conv<ZZ>(6666);
   aa[7] = conv<ZZ>(7777);
   aa[8] = conv<ZZ>(9272);
   testMWCProposal<Int, Real>("Marsa 16b", 8, 16, aa, maxdim);
   // We confirm here that m and (m-1)/2 are not prime.
   Int m = conv<ZZ>("3155138487111751905571868744270142781194239");
   IntFactorization<Int> fact(m);
   fact.factorize();
   std::cout << "Factorization of m = " << fact.toString() << "\n";
   fact.setNumber(m - Int(1));
   fact.factorize();
   std::cout << "Factorization of m - 1 = " << fact.toString() << "\n";
   // For this one,
   // m   = 3155138487111751905571868744270142781194239
   //     = 517854180589 * 6092716068301586638428281517851
   // m-1 = 3155138487111751905571868744270142781194238
   //     = 2 * 107385419 * 14690721126262737334813904037894301

   // Example 3 from Marsaglia (1994).
   aa.SetLength(3);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>("1111111464");
   aa[2] = conv<ZZ>("1111111464");
   // fomdual.setVerbosity(4);
   testMWCProposal<Int, Real>("Marsa 32", 2, 32, aa, maxdim);

   // Int b = getPow2<Int>(32);
   //m = computeLCGModulusMWC(b, aa);
   //std::cout << "11404 * 1111111464 % m = " << Int(11407) * aa[1] % m << "\n";

   // ==============================================================================

   std::cout << "\n=============================================================\n";
   std::cout << "MWC generator examined by Couture and L'Ecuyer (1997).\n\n";
   maxdim = 16;
   // fomdual.setVerbosity(4, 0.1);

   // Example from Couture and L'Ecuyer (1997).
   aa.SetLength(9);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>(14);
   aa[2] = conv<ZZ>(18);
   aa[3] = conv<ZZ>(144);
   aa[4] = conv<ZZ>(1499);
   aa[5] = conv<ZZ>(2083);
   aa[6] = conv<ZZ>(5273);
   aa[7] = conv<ZZ>(10550);
   aa[8] = conv<ZZ>(45539);
   testMWCProposal<Int, Real>("Cout-Lec 16", 8, 16, aa, maxdim);

   // ==============================================================================

   std::cout << "\n=============================================================\n";
   std::cout << "Some MWC generator proposed by Goresky and Klapper (2003).\n\n";

   // b = 2^{24}, k = 48.
   // Here we temporarily use Roger's bounds because we need more than 48 dimensions.
   maxdim = 52;
   NormaRogers normaDualR(52);
   fomdual.setNormalizer(normaDualR);

   aa.SetLength(49);
   clear(aa);
   aa[0] = 3;
   aa[14] = -2;
   aa[38] = -2;
   aa[46] = -2;
   aa[48] = 2;
   testMWCProposal<Int, Real>("GK 24a", 48, 24, aa, maxdim);

   // b = 2^{24}, k = 41.
   maxdim = 48;
   fomdual.setNormalizer(normaDual);
   aa.SetLength(42);
   clear(aa);
   aa[0] = 3;
   aa[14] = -4;
   aa[38] = -2;
   aa[41] = 2;
   testMWCProposal<Int, Real>("GK 24b", 41, 24, aa, maxdim);

   // b = 2^{25}, k = 22.
   aa.SetLength(23);
   clear(aa);
   aa[0] = 3;
   aa[4] = 2;
   aa[6] = -2;
   aa[11] = 2;
   aa[15] = 2;
   aa[16] = -2;
   aa[17] = -2;
   aa[20] = -2;
   aa[22] = 2;
   testMWCProposal<Int, Real>("GK 25", 22, 25, aa, maxdim);

   // b = 2^{32}, k = 33.
   aa.SetLength(34);
   clear(aa);
   aa[0] = 5;
   aa[4] = -4;
   aa[11] = -4;
   aa[14] = -4;
   aa[20] = -4;
   aa[33] = 4;
   testMWCProposal<Int, Real>("GK 32", 33, 32, aa, maxdim);

//   return 0;

   // ==============================================================================

   std::cout << "\n=============================================================\n";
   std::cout << "MWC generator proposed by Vigna (2021), on the Internet.\n\n";
   maxdim = 8;

   // MWC128.c  $k=1$, $a_1 = \texttt{0xffebb71d94fcdaf9} = 18441034436880161529$.
   aa.SetLength(2);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>(0xffebb71d94fcdaf9);
   // aa[1] = conv<ZZ>("18441034436880161529");
   testMWCProposal<Int, Real>("Vigna MWC128", 1, 64, aa, maxdim);

   // MWC192.c  $k=2$, $a_2 = {0xffa04e67b3c95d86} = 18419808683250244998$, $a_1=0$.
   aa.SetLength(3);
   aa[1] = Int(0);
   aa[2] = conv<ZZ>(0xffa04e67b3c95d86);
   testMWCProposal<Int, Real>("Vigna MWC192", 2, 64, aa, maxdim);

   Coordinates coord( { 1, 3, 4 });
   IntLattice<Int, Real> proj(lcg.getModulus(), 4); // Lattice used for projections.
   lcg.buildProjectionDual(proj, coord);
   proj.dualize();
   std::cout << "Spectral test for the worst projection in dual lattice:\n";
   fomdual.computeMeritOneProj(proj, coord);
   proj.dualize();

   // MWC256.c  $k=3$, $a_3 = {0xfff62cf2ccc0cdaf}$, $a_1=a_2=0$.
   aa.SetLength(4);
   aa[1] = 0;
   aa[2] = 0;
   aa[3] = conv<ZZ>(0xfff62cf2ccc0cdaf);
   testMWCProposal<Int, Real>("Vigna MWC256", 3, 64, aa, maxdim);

   coord.erase(3);
   coord.insert(5);
   proj.setModulus(lcg.getModulus());
   lcg.buildProjectionDual(proj, coord);
   proj.dualize();
   std::cout << "Spectral test for the worst projection in dual lattice:\n";
   fomdual.computeMeritOneProj(proj, coord);
   proj.dualize();

   // GMWC128.c  $k=1$, $a_1 = {0xff002aae7d81a646} = 18374733408589948486$,
   //              $a_0^{-1} = {0x9b1eea3792a42c61} = 11177628849584483425$,
   // m = 0xff002aae7d81a646007d084a4d80885f
   //   = 338954004610899541305165203194907756639
   aa.SetLength(2);
   // aa[0] = conv<ZZ>("35193487309703263");
   // aa[1] = conv<ZZ>("18374733408589948486");
   aa[0] = conv<ZZ>(0x7d084a4d80885f);
   aa[1] = conv<ZZ>(0xff002aae7d81a646);
   testMWCProposal<Int, Real>("Vigna GMWC128", 1, 64, aa, maxdim);

   // GMWC256.c  $k=3$, $a_3 = 0xff963a86efd088a2$,
   //  $a_0 = 0x54c3da46afb70f$, $a_0^{-1} = 0xbbf397e9a69da811.
   // m = 0xff963a86efd088a2000000000000000000000000000000000054c3da46afb70f
   //   = 115605207387626077441750929865143416680831684312500419426686549709778808977167
   aa.SetLength(4);
   aa[0] = conv<ZZ>(0x54c3da46afb70f);
   aa[1] = 0;
   aa[2] = 0;
   aa[3] = conv<ZZ>(0xff963a86efd088a2);
   testMWCProposal<Int, Real>("Vigna GMWC256", 3, 64, aa, maxdim);

   coord.insert(2);
   proj.setModulus(lcg.getModulus());
   lcg.buildProjectionDual(proj, coord);
   proj.dualize();
   std::cout << "Spectral test for the worst projection in dual lattice:\n";
   fomdual.computeMeritOneProj(proj, coord);
   proj.dualize();

   return 0;

   // ==============================================================================

   // Jump ahead.

   int64_t k = 1;
   int64_t e = 64;
   aa.SetLength(k + 1);
   aa[0] = Int(-1);
   aa[1] = conv<ZZ>(0xffebb71d94fcdaf9);
   Int b = NTL::power(Int(2), e);
   m = computeLCGModulusMWC(b, aa);
   IntVec xx;
   xx.SetLength(k);
   Int c;
   Int s0;

   std::cout << "____________________________________________________\n";
   std::cout << "MWC to LCG to MWC \n\n";
   std::cout << "MWC of order k = " << k << ",  modulo b = 2^" << e << " = " << b << "\n";
   std::cout << "Coefficients aa = " << aa << "\n";

   c = 12345;
   xx[0] = 123456789;
   MWCtoLCGState<Int>(s0, b, aa, xx, c);
   std::cout << "xx = " << xx << "\n";
   std::cout << "c = " << c << "\n";
   std::cout << "State s0 = " << s0 << "\n\n";

   // s0 = conv<Int>("1234567890");
   LCGtoMWCState<Int>(xx, c, b, aa, s0);
   std::cout << "State s0 = " << s0 << "\n";
   std::cout << "xx = " << xx << "\n";
   std::cout << "c = " << c << "\n";
   MWCtoLCGState<Int>(s0, b, aa, xx, c);
   std::cout << "State s0 = " << s0 << "\n";

   return 0;
}

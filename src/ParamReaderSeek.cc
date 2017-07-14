/**
ParamReaderSeek.cc for ISO C++
version
modified
authors: Hicham Wahbi
         Frederik Rozon
         Richard Simardr
*/

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latmrg/ParamReaderSeek.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>



using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

//===========================================================================

ParamReaderSeek::ParamReaderSeek ()
{}


//===========================================================================

ParamReaderSeek::ParamReaderSeek (string fileName): ParamReader (fileName)
{}


//===========================================================================

ParamReaderSeek::~ParamReaderSeek ()
{
}


//===========================================================================

void ParamReaderSeek::read (SeekConfig & config)
{
   getLines ();
   int ln = 0;
   unsigned int lnu;
   int s;

   readCriterionType (config.criter, ++ln, 0);
   if (config.criter == SPECTRAL)
      readNormaType (config.normaType, ln, 1);

   readBool (config.readGenFile, ++ln, 0);
   if (config.readGenFile)
      readString (config.fileName, ln, 1);
   readInt (config.J, ++ln, 0);

   config.Ms = new MScal[config.J];
   config.Ks = new int[config.J];

   config.compon.reserve (config.J);
   int maxOrder = 0;
   Component comp;
   for (int i = 0; i < config.J; i++) {
      readGenType (comp.genType, ++ln, 0);
      readNumber3 (comp.modulus.m, comp.modulus.b, comp.modulus.e,
                   comp.modulus.c, ++ln, 0);
      if (0 == comp.modulus.e)
         comp.modulus.init (comp.modulus.m);
      else
         comp.modulus.init (comp.modulus.b, comp.modulus.e, comp.modulus.c);
      config.Ms[i] = comp.modulus.m;

      readInt (comp.k, ++ln, 0);
      if (maxOrder < comp.k)
         maxOrder = comp.k;
      config.Ks[i] = comp.k;
      readBool (comp.PerMax, ++ln, 0);

      readDecompType (comp.F1, ++ln, 0);
      if (comp.F1 == DECOMP_WRITE || comp.F1 == DECOMP_READ) {
         comp.file1.reserve (MAX_WORD_SIZE);
         readString (comp.file1, ln, 1);
      }
      readDecompType (comp.F2, ++ln, 0);
      if (comp.F2 == DECOMP_WRITE || comp.F2 == DECOMP_READ) {
         comp.file2.reserve (MAX_WORD_SIZE);
         readString (comp.file2, ln, 1);
      }

      readImplemCond (comp.implemCond, ++ln, 0);
      if (POWER_TWO == comp.implemCond) {
         readInt (comp.NumBits, ln, 1);
         readInt (comp.HighestBit, ln, 2);
      } else if (EQUAL_COEF == comp.implemCond ||
                 ZERO_COEF == comp.implemCond) {
         readInt (comp.ncoef, ln, 1);
         comp.Icoef = new int[1 + comp.ncoef];
         for (s = 1; s < comp.ncoef; s++)
            readInt (comp.Icoef[s], ln, 1 + s);
         comp.Icoef[0] = 0;
         comp.Icoef[comp.ncoef] = comp.k;
      }

      comp.b.SetLength (1 + comp.k);
      comp.c.SetLength (1 + comp.k);

      lnu = ln;
      if (ZERO_COEF == comp.implemCond) {
         SetZero (comp.b, comp.k);
         SetZero (comp.c, comp.k);
         readInterval (comp.b, comp.c, lnu, comp.ncoef);
         for (s = comp.ncoef; s >= 1; s--) {
            comp.b[comp.Icoef[s]] = comp.b[s];
            comp.c[comp.Icoef[s]] = comp.c[s];
            if (s < comp.Icoef[s]) {
               comp.b[s] = 0;
               comp.c[s] = 0;
            }
         }

      } else if (EQUAL_COEF == comp.implemCond) {
         readInterval (comp.b, comp.c, lnu, comp.ncoef);
         for (s = comp.ncoef; s >= 1; s--) {
            for (int i = comp.Icoef[s]; i > comp.Icoef[s-1]; i--) {
               comp.b[i] = comp.b[s];
               comp.c[i] = comp.c[s];
            }
         }
      } else
         readInterval (comp.b, comp.c, lnu, comp.k);
      ln = lnu;

      if (comp.genType == RANK1)
         comp.b[1] = comp.c[1] = 1;
      checkBound (comp.modulus.m, comp.b, comp.k);
      checkBound (comp.modulus.m, comp.c, comp.k);

      if (POWER_TWO == comp.implemCond) {
         for (int i = 1; i <= comp.k; i++) {
            comp.b[i] = 0;
            comp.c[i] = comp.HighestBit;
         }
      }
      readSearchMethod (comp.searchMethod, ++ln, 0);
      if (RANDOM == comp.searchMethod) {
         readInt (comp.numReg, ln, 1);

         readInt (comp.H, ln, 2);
         readInt (comp.Hk, ln, 3);
      } else
         comp.numReg = 1;

      config.compon.push_back (comp);
   }

   readInt (config.C, ++ln, 0);
   readDoubleVect (config.minMerit, ++ln, 0, config.C, 0);
   readDoubleVect (config.maxMerit, ++ln, 0, config.C, 0);
   readIntVect (config.numGen, ++ln, 0, config.C, 0);
   readInt (config.d, ++ln, 0);
   if (config.d < 1)
      MyExit (1, "ParamReaderSeek:   config.d < 1");
   config.td = new int[1 + config.d];
   readIntVect (config.td, ++ln, 0, 1 + config.d, 0);
   if (config.td[0] <= maxOrder)
      config.td[0] = maxOrder + 1;

   readBool (config.dualF, ++ln, 0);
   readLatticeType (config.latType, ++ln, 0);

   readInt (config.lacGroupSize, ++ln, 0);
   readInt (config.lacSpacing, ln, 1);
   readLong (config.maxNodesBB, ++ln, 0);

   double tem;
   readDouble (tem, ++ln, 0);
   char c;
   readChar (c, ln, 1);
   switch (c) {
   case 's':
      config.duration = tem;
      break;
   case 'm':
      config.duration = tem * 60;
      break;
   case 'h':
      config.duration = tem * 3600;
      break;
   case 'd':
      config.duration = tem * 86400;
      break;
   default:
      config.duration = 0;
      cout << "***** reading duration:  IMPOSSIBLE CASE";
   }

   readLong (config.seed, ++ln, 0); // seed of the generator
//   readLong (config.s2, ln, 2);
   readOutputType (config.outputType, ++ln, 0);
}

}

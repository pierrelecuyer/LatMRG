#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "latticetester/Rank1Lattice.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaMinkowski.h"
#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/WriterRes.h"
#include "latticetester/Types.h"

#include "latmrg/LatTestBeyer.h"
#include "latmrg/LatTestSpectral.h"
#include "latmrg/LatTestPalpha.h"
#include "latmrg/ReportHeader.h"
#include "latmrg/Merit.h"
#include "latmrg/TestProjections.h"
#include "latmrg/Zone.h"
#include "latmrg/Chrono.h"
#include "latmrg/Const.h"
#include "latmrg/SeekConfig.h"
#include "latmrg/ParamReaderSeek.h"
#include "latmrg/IntPrimitivity.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/KorobovLattice.h"
#include "latmrg/MRGLatticeFactory.h"
#include "latmrg/MRGComponent.h"

using namespace std;
// using namespace NTL;
using namespace LatMRG;
using namespace LatticeTester;


namespace
{
  //===========================================================================

  Normalizer<NScal> *normal;
  // Flag to ensure that Normalizer are created only once;
  // must know m and k for that.
  bool isFirstTest = true;
  NormType Norm = L2NORM;

  SeekConfig<MScal> config;
  ofstream fout;
  string outfile;
  // streambuf *psbuf;
  Chrono timer;

  MRGComponent<MScal> **compJ;
  MMat coef;
  LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *lattice;
  LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *master;
  LatticeTest<MScal, NScal> *latTest;
  LatticeTest<MScal, NScal> **pool;            // vecteur de LatticeTest*
  int poolLen = 0;               // longueur du vecteur pool <= numGen

  // What are those? Nobody knows.
  Zone<MScal> **zone;
  Zone<MScal> **reg;
  IntPrimitivity<MScal> **primJ;

  long numATried = 0;
  // long numAPrimitive = 0;
  long *CoNoElemPrim;
  long *CoElemPrim;
  long *TotEP;
  long *CoPerMax;
  long *CoRegions;
  long *Hk;
  long *H;
  double *minVal;
  // long ExceedBB = 0;    // Number of cases for which m_countNodes > maxNodesBB

  void SeekGen (int);



  //===========================================================================

  void SortBestGen ()
  {
    for (int r = 0; r < poolLen - 1; r++) {
      LatticeTest<MScal, NScal> *latTest1 = pool[r];
      double merit1 = latTest1->getMerit ().getWorstMerit ();
      for (int s = r + 1; s < poolLen; s++) {
        LatticeTest<MScal, NScal> *latTest2 = pool[s];
        double merit2 = latTest2->getMerit ().getWorstMerit ();
        if (merit2 > merit1) {
          latTest = pool[s];
          pool[s] = pool[r];
          pool[r] = latTest;
          merit1 = merit2;
        }
      }
    }
  }


  //===========================================================================

  int readConfigFile (int argc, char **argv)
  {
    if (argc != 2) {
      cerr << "\n***  Usage: " << argv[0] << " <data file>" << endl;
      return -1;
    }
    // Lecture des paramètres
    string fname (argv[1]);
    //string fname ("./inputTestFiles/seekZZDD_test1");
    fname += ".dat";
    ParamReaderSeek<MScal, NScal> paramRdr (fname.c_str ());
    paramRdr.read (config);

    // Writer *rw;
    switch (config.outputType) {
      case RES:
        fname = argv[1];
        //fname = "./inputTestFiles/seekZZDD_test1";
        fname += ".res";
        // rw = new WriterRes (fname.c_str());
        outfile = fname;
        break;
      case TEX:
        fname = argv[1];
        fname += ".tex";
        return -2;
        // rw = new WriterTex(fname.c_str());
        break;
      case TERM:
        cerr << "\n*** SeekMain::readConfigFile::outputType: " <<
          "TERM is not implemented." << endl;
        throw std::invalid_argument ("SeekMain::readConfigFile");
        // rw = new WriterRes (&cout);
        // psbuf = cout.rdbuf(); // get cout streambuf
        // fout.rdbuf(psbuf); // assign streambuf to fout
        break;
      case GEN:
        fname = argv[1];
        fname += ".gen";
        outfile = fname;
        break;
      default:
        cerr << "\n*** outputType:   no such case" << endl;
        return -2;
    }
    fname.clear ();
    return 0;
  }


  //===========================================================================

  /*
   *void debug (int dim1, int dim2)
   *{
   *   bool racF = false;
   *   bool invertF = true;
   *   cout <<
   *   "=================================================================\n\n";
   *   for (int j = 0; j < poolLen; j++) {
   *      LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *lat = pool[j]->getLattice ();
   *      cout << "A = "; lat->write();
   *      int dimWorst;
   *      double S_T = pool[j]->getMerit ().getST (dim1, dim2, dimWorst);
   *      if (lat->getNorm () == L2NORM)
   *         racF = true;
   *      cout << "\t\tS_" << dimWorst << " = " << S_T << endl << endl;
   *      cout << "   dim\t\t\tl_dim\t\t     S_dim" << endl;
   *      cout << pool[j]->getMerit ().toString (dim1, dim2, racF,
   *                                             invertF) << endl;
   *   }
   *}
   *
   */

  //===========================================================================
  /*
     void printPool (char *mess)
     {
     fout << "--------------------------------------------------" <<
     mess << endl;
     for (int i = 0; i < poolLen; i++) {
     LatticeTest *latTest = pool[i];
     LatticeTester::IntLattice<MScal, BScal, BVect, BMat, NScal, NVect, RScal> *lat = latTest->getLattice ();
     int dimWorst;
     double S_T = latTest->getMerit ().getST (config.td[0],
     config.td[1], dimWorst);
     if (lat->getNorm () == L2NORM) {
     S_T = sqrt (S_T);
     }
     if (config.J > 1) {
     for (int j = 0; j < config.J; j++) {
     fout << " " << toString (lat->comp[j]->a, 1, lat->comp[j]->k);
     if (j < config.J - 1)
     fout << endl;
     }
     } else {
     fout << " "; lat->write();
     }
     fout << "\t\tS_" << dimWorst << " = " << S_T << endl;
     }
     }

*/
  //===========================================================================

  void PrintComponentData (const Component<MScal> & comp, int j)
  {
    // Prints parameter of component j of combined generator
    fout << "Component " << j + 1 << endl;
    fout << "-----------\n";
    fout << "   GenType                : " << toStringGen (comp.
        genType) << endl;
    fout << "   Modulus m              : " << comp.modulus.m << endl;
    fout << "   Order k                : " << comp.k << endl;
    // fout << " n_j : " << comp.n << endl;
    // fout << " rho_j = (m_j)^k_j-1 : " << comp.rho << endl;
    if (comp.PerMax) {
      fout << "   Factors of m-1         : " <<
        compJ[j]->ifm1.toString ();
      fout << "   Factors of r           : " <<
        compJ[j]->ifr.toString ();
    }
    fout << "   Search method                      : " <<
      toStringSearchMethod (comp.searchMethod) << endl;
    string margin = "   Bounds : ";
    for (int i = 0; i < comp.k; i++) {
      if (i != 0)
        margin = "            ";
      fout << margin << "a" << i+1 << " from : " << comp.b[i] << endl;
      fout << "                 to : " << comp.c[i] << endl;
    }
    fout << "   Implementation condition           : " <<
      toStringImplemCond (comp.implemCond) << endl;
    fout << "   Maximum period required            : ";
    if (comp.PerMax)
      fout << "yes";
    else
      fout << "no";
    fout << endl << "\n-----------------------------------------------" <<
      endl;
  }


  //===========================================================================

  void PrintPolyStats (const Component<MScal> & comp, int j)
  {
    fout << "Component " << j + 1 << endl;
    fout << "-----------\n";
    if (comp.PerMax) {
      fout << "    Values of a_" << comp.k << " tried                   : "
        << CoNoElemPrim[j] + CoElemPrim[j] << endl;
      fout << "    Values of a_" << comp.k <<
        " primitive element       : " << CoElemPrim[j] << endl;
      fout << "    Polynomials with a_" << comp.k <<
        " primitive element: " << TotEP[j] << endl;
      fout << "    Primitive polynomials                 : "
        << CoPerMax[j] << endl;
    } else
      fout << "    Num. of polynomials examined          : " <<
        TotEP[j] << endl;
    fout << "\n-----------------------------------------------" << endl;
  }


  //===========================================================================

  void PrintResults ()
  {
    fout.open (outfile.c_str());
    Component<MScal> & comp0 = config.compon[0];

    fout << "SEARCH for good ";
    if (comp0.genType == KOROBOV)
      fout << "KOROBOVs" << endl;
    else if (comp0.genType == RANK1)
      fout << "RANK1 Lattices of order " << comp0.k << endl;
    else
      fout << "MRGs of order " << comp0.k << endl;

    // fout << "-----------------------------------------------\n";
    fout << "\nDATA \n";
    fout << "-----------------------------------------------\n";
    for (int j = 0; j < config.J; j++)
      PrintComponentData (config.compon[j], j);

    fout << "   Test                               : "
      << toStringCriterion (config.criter) << endl;
    fout << "   Normalizer                         : "
      << toStringNorma (config.normaType) << endl;
    fout << "   num categories C                   : " << config.C << endl;
    fout << "   min merit                          : " <<
      toString < double *>(config.minMerit, 0, config.C) << endl;
    fout << "   max merit                          : " <<
      toString < double *>(config.maxMerit, 0, config.C) << endl;
    fout << "   d                                  : " << config.d << endl;
    fout << "   td                                 : " <<
      toString (config.td, 0, config.d+1) << endl;

    if (RANDOM == comp0.searchMethod) {
      fout << "   Seed for RNG                       : " <<
        config.seed << endl;
    }
    fout << "   Max nodes in branch-and-bound      : " <<
      config.maxNodesBB << endl;
    fout << "   Lattice Type                       : " <<
      toStringLattice (config.latType) << endl;
    fout << "   Lattice                            : ";
    if (config.dualF)
      fout << "DUAL" << endl;
    else
      fout << "PRIMAL" << endl;
    fout << "   Speed option                       : " << config.speed << endl;
    fout << "\n\nRESULTS" << endl;
    fout << "-----------------------------------------------\n";
    for (int j = 0; j < -config.J; j++)
      PrintPolyStats (config.compon[j], j);

    fout << "Num. of generators tested             : " <<
      numATried << endl;
    fout << "Num. of generators conserved          : " <<
      poolLen << endl;
    fout << "Total CPU time (after setup)          : " <<
      timer.toString () << endl;
    fout << setprecision (5) << endl << endl;

    bool rac = false;
    if (Norm == L2NORM)
      rac = true;

    for (int i = 0; i < config.C; i++) {
      fout <<
        "+--------------------------------------------------------------------"
        << endl;
      fout << "   " << poolLen << " generators retained for criterion M_";
      for (int j = 1; j < config.d; j++)
        fout << config.td[j] << ",";
      fout << config.td[config.d] << endl;
      fout <<
        "+--------------------------------------------------------------------"
        << endl;
      fout << "   A\t\t\tM_*" << endl;
      fout <<
        "+--------------------------------------------------------------------"
        << endl;
      for (int s = 0; s < poolLen; s++) {
        LatticeTest<MScal, NScal> *latTest = pool[s];
        LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *lat = latTest->getLattice ();
        int dimWorst = latTest->getMerit ().getDimWorst ();
        double S_T = latTest->getMerit ().getWorstMerit ();
        if (rac)
          S_T = sqrt (S_T);
        if (config.J > 1) {
          fout << "Cannot print components because things have been done \
            badly\n";
          /*
             for (int j = 0; j < config.J; j++) {
             fout << " " << toString (lat->comp[j]->a, 1,
             lat->comp[j]->k);
             if (j < config.J - 1)
             fout << endl;
             }
             */
        } else {
          fout << " [" << lat->toStringCoef ();
        }
        fout << "]\n S_";
        if (1 == config.d)
          fout << dimWorst;
        else
          fout << "*";
        fout << " = " << S_T << "\n" << endl;
        if (1 == config.d) {
          fout << "   dim\t\t\tl_dim\t\t     S_dim" << endl;
          fout << latTest->getMerit ().toString (config.td[0],
              config.td[1], rac, config.invertF);
        }
        /* if (config.dualF) fout << lat->getDualBasis ().toString(); else
           fout << "m*" << lat->getPrimalBasis ().toString(); */
        fout <<
          "\n+-----------------------------------------------------------------"
          << endl;
      }
    }
    fout.close();
  }

  //===========================================================================

  // The output function for the outputType GEN. This function prints the 
  // coefficients of the retained generators in a .gen file. One generator is
  // printed per line, coefficients are separeted by a blank caracter. No 
  // figure of merit is printed.
  //
  // This currently do not support J > 1
  void PrintGen () {
    fout.open (outfile.c_str());
    MVect coefs;
    fout << poolLen << std::endl;
    for (int i = 0; i < config.C; i++) {
      for (int s = 0; s < poolLen; s++) {
        LatticeTest<MScal, NScal> *latTest = pool[s];
        LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *lat = latTest->getLattice ();
        if (config.J > 1) {
          fout << "Cannot print components because things have been done \
            badly\n";
          /*
             for (int j = 0; j < config.J; j++) {
             fout << " " << toString (lat->comp[j]->a, 1,
             lat->comp[j]->k);
             if (j < config.J - 1)
             fout << endl;
             }
             */
        } else {
          fout << lat->toStringCoef () << std::endl;
        }
      }
    }
    fout.close();
  }

  //===========================================================================

  // This will perform the tests specified by the user in the .dat file.
  void Test ()
  {
  bool stationary = true;
  Component<MScal> & comp0 = config.compon[0];

  if (config.J > 1) {         // On doit faire la combinaison des MRG
    for (int s = 0; s < config.J; s++)
      compJ[s]->setA (coef[s]);
    lattice = MRGLatticeFactory<MScal, NScal>::fromCombMRG (compJ, config.J,
        config.getMaxDim (), 0, config.latType, Norm);

  } else {
    if (comp0.genType == MRG || comp0.genType == LCG)
      lattice =
        new MRGLattice<MScal, NScal> (comp0.modulus.m, coef[0],
            config.getMaxDim (), comp0.k, config.latType, Norm);
    else if (comp0.genType == KOROBOV)
      lattice =
        new KorobovLattice<MScal, NScal> (comp0.modulus.m, coef[0][0],
            config.getMaxDim (), config.latType, Norm);
    else if (comp0.genType == RANK1) {
      stationary = false;
      lattice = new LatticeTester::Rank1Lattice<MScal, BScal, NScal, RScal> (
          comp0.modulus.m, coef[0], config.getMaxDim (), Norm);
    }
  }

  if (isFirstTest) {
    normal = lattice->getNormalizer (config.normaType, config.alpha, config.dualF);
    isFirstTest = false;
  }

  lattice->buildBasis (config.td[0]);

  switch (config.criter) {
    case SPECTRAL:
      latTest = new LatTestSpectral<MScal, NScal> (normal, lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setInvertFlag (config.invertF);
      latTest->setMaxAllDimFlag (true);
      if (1 == config.d) {
        if (config.speed == 0) {
          latTest->test (config.td[0], config.td[1], minVal);
          latTest->getMerit ().getST (config.td[0], config.td[1]);
        } else if (config.speed == 1) {
          latTest->quicktest (config.td[0], config.td[1], minVal, 1);
          latTest->getMerit ().getST (config.td[0], config.td[1]);
        } else if (config.speed == 2) {
          latTest->quicktest (config.td[0], config.td[1], minVal, 2);
          latTest->getMerit ().getST (config.td[0], config.td[1]);
        }
      } else {
        if (comp0.genType == MRG || comp0.genType == LCG)
          master = new MRGLattice<MScal, NScal> (*(MRGLattice<MScal, NScal> *) lattice);
        else if (comp0.genType == KOROBOV)
          master = new KorobovLattice<MScal, NScal> (*(KorobovLattice<MScal, NScal> *) lattice);
        else if (comp0.genType == RANK1)
          master = new LatticeTester::Rank1Lattice<MScal, BScal, NScal, RScal>
            (*(LatticeTester::Rank1Lattice<MScal, BScal, NScal, RScal> *) lattice);

        master->buildBasis (config.td[1]);
        TestProjections<MScal, NScal> proj (master, lattice, latTest, config.td, config.d);
        proj.setDualFlag (config.dualF);
        //proj.setPrintF (true);
        double merit = proj.run (stationary, false, minVal);
        latTest->getMerit ().setWorstMerit (merit);
        delete master;
      }
      break;

    case BEYER:
      latTest = new LatTestBeyer<MScal, NScal>  (lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setInvertFlag (config.invertF);
      latTest->setMaxAllDimFlag (true);
      latTest->test (config.td[0], config.td[1], minVal);
      break;

    case PALPHA:
      latTest = new LatTestPalpha<MScal, NScal> (normal, lattice);
      latTest->setDualFlag (config.dualF);
      latTest->setMaxAllDimFlag (true);
      latTest->test (config.td[0], config.td[1], minVal);
      break;

    default:
      fout << "Invalid CriterionType " << config.criter << endl;
      exit (1);
      break;
  }

}


//===========================================================================
// I don't know what purpose this had.
#if 0

void VerifyCategories (BVect & Me2)
{
  // Called by InsideExam.  Me2 is either Q2 or S2.
  int MaxI;
  double ToBeat;
  int Ind = Order + 1;
  int IndMin = Ind;
  for (int C = 0; C < NbCat; ++C) {
    ++NumTried[C];
    if ((C + 1 < NbCat) || ((C + 1 == NbCat) && TestCompleted)) {
      // On trouve ou le min est atteint pour cette categorie.
      if ((Criterion == Spectral) && (TabDim[C] <= MaxDimS2))
        MaxI = TabDim[C];
      else
        MaxI = MaxDimS2;
      while (Ind <= MaxI) {
        if (Me2[Ind] <= 0.0)
          return ;
        if (Me2[Ind] < Me2[IndMin])
          IndMin = Ind;
        ++Ind;
      }
      if (Me2[IndMin] > MaxMerit[C] * MaxMerit[C])
        return ;
      WITH PireGen[C] ^ DO
        if ((DimMerit == 0) || ((MaximCat[C]) && (Me2[IndMin] >= Merit))
            || ((!MaximCat[C]) && (Me2[IndMin] <= Merit)))
          ConserverGen (C, IndMin);

    }
  }
}
#endif

//===========================================================================

/* Find the minimum merit amongst the gen in category 0, i.e. in pool.
 * Returns this merit.
 * Side effect: swap that gen with the gen in element 0, so that on return,
 * the worst gen is in pool[0] */
double FindMinMerit ()
{
  double minMerit = 1.0e300;
  int minj = -1;
  for (int j = 0; j < poolLen; j++) {
    double curMerit = pool[j]->getMerit ().getWorstMerit ();
    if (curMerit < minMerit) {
      minj = j;
      minMerit = curMerit;
    }
  }

  if (minj > 0) {
    LatticeTest<MScal, NScal> *t = pool[minj];
    pool[minj] = pool[0];
    pool[0] = t;
  }

  return minMerit;
}


//===========================================================================

// It think tries to fit a new generator that succeeded in the tests in the 
// table of the best generators we found and orders it.
void CompareMerit ()
{
  int s;
  if (poolLen < config.numGen[0]) {
    pool[poolLen] = latTest;
    poolLen++;
    FindMinMerit ();

  } else {
    double minMerit = pool[0]->getMerit ().getWorstMerit ();
    double curMerit = latTest->getMerit ().getWorstMerit ();

    if (minMerit < curMerit) {
      LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *lat = pool[0]->getLattice ();
      delete lat;
      delete pool[0];
      pool[0] = latTest;
      lat = latTest->getLattice ();
      minMerit = FindMinMerit ();
      int dim1 = config.td[0];
      int dim2 = config.td[1];
      for (s = dim1; s <= dim2; s++)
        minVal[s] = minMerit;

    } else {
      LatticeTester::IntLattice<MScal, BScal, NScal, RScal> *lat = latTest->getLattice ();
      delete lat;
      delete latTest;
    }
  }
}


//===========================================================================

// This should be called ExamThisai but w/e. This examines the component a_i
// (that was just set by the function that called this one) of the j-th MRG.
// If 1<i<k this will call a function that will increment a_{i-1}. If i=1, if
// j<J we call a function that will examine the next MRG component. If j=J, we
// call a function that will perform the spectral test with the set of 
// components currently in place.
void ExamThisaj (int j, int i, bool Pow2, ProcII Exam)
{
  // Called by InsideExam.  Examines the current a_j
  Component<MScal> & comp = config.compon[j];
  if ((i == comp.k) && (coef[j][i] == 0))
    return ;

  bool Ok = true;
  if (!Pow2 && comp.PerMax && (i == comp.k)) {
    // On verifie cond.(1) pour periode max
    if (comp.modulus.primeF) { // m_j is a prime
      if (primJ[j]->isPrimitiveElement (coef[j], comp.k))
        ++CoElemPrim[j];
      else {
        Ok = false;
        ++CoNoElemPrim[j];
      }
    } else {
      Ok = comp.modulus.perMaxPowPrime (coef[j][i]);
      if (Ok)
        ++CoElemPrim[j];
      else
        ++CoNoElemPrim[j];
    }
  }

  if (!Ok)
    return ;

  if (comp.implemCond == EQUAL_COEF) {
    int r = 1;
    while ((r < comp.ncoef) && (i > comp.Icoef[r]))
      ++r;
    r--;
    int s = i - 1;
    while (s > comp.Icoef[r]) {
      coef[j][s] = coef[j][i];
      s--;
    }
    i = s + 1;

  } else if (comp.implemCond == ZERO_COEF) {
    int r = comp.ncoef;
    while ((r > 0) && (i <= comp.Icoef[r]))
      --r;
    i = comp.Icoef[r] + 1;
    if (i < 1)
      i = 1;
  }

  if (i == 0) {
    ++TotEP[j];
    if (Pow2 || (!comp.PerMax) || (comp.k == 1) ||
        compJ[j]->maxPeriod23 (coef[j])) {
      // Note: MaxPeriod verifie cond. (2 et 3) pour per. max.
      // If power of prime, then j = i = Orderj [j] = 1.
      ++CoPerMax[j];
      SeekGen (j + 1);
    }

  } else {
    Exam (j, i - 1);
  }
}


//===========================================================================

// ??
void ExamBits (MScal q, int j, int i, int b0, int b1, int NbBits,
    bool Pow2, ProcII Exam)
{
  /* Called by InsideExam (and calls itself recursively) to examine all
     bit patterns with less than NbMaxBits [j] in total, that can be
     obtained by switching from 0 to 1 some of the bits of q from bit b0
     to bit b1. */
  Component<MScal> & comp = config.compon[j];

  for (int b = b0; b <= b1; b++) {
    coef[j][i] = q + TWO_EXP[b];
    ExamThisaj (j, i, Pow2, Exam);
    if ((b > b0) && (NbBits < comp.NumBits))
      ExamBits (coef[j][i], j, i, b0, b - 1, NbBits + 1, Pow2, Exam);
    coef[j][i] = coef[j][i] - TWO_EXP[b];
    coef[j][i] = coef[j][i] - TWO_EXP[b];
    ExamThisaj (j, i, Pow2, Exam);
    if ((b > b0) && (NbBits < comp.NumBits))
      ExamBits (coef[j][i], j, i, b0, b - 1, NbBits + 1, Pow2, Exam);
    coef[j][i] = coef[j][i] + TWO_EXP[b];
  }
}


//===========================================================================

// Increments the coefficient that Z refers to. At each step this calls a 
// function that will look at the next coefficient, or that will call a method
// to test the current state of the generator.
void InsideExam (Zone<MScal> * Z, int j, int i, ProcII exam)
{
  /*
   *  Examines all the possible values for coefficient $a_i$ of component
   *  `j` in this region according to the chosen criteria, and
   *  calls different methods depending on the criteria. \texttt{exam} is
   *  the (recursive) method that called \texttt{InsideExam}. */
  Component<MScal> & comp = config.compon[j];
  MScal q, Temp;
  bool Pow2 = false;
  MScal Eight;
  Eight = 8;
  q = Z->getInf ();

  if (comp.PerMax && (!comp.modulus.primeF) && (comp.modulus.b == 2)
      && comp.modulus.c == 0) {
    Pow2 = true;             // m is a power of 2: want a mod 8 = 5 only.
    Modulo (q, Eight, Temp);
    while (Temp != 5) {
      ++q;
      Modulo (q, Eight, Temp);
    }
  }

  if (comp.implemCond == POWER_TWO) {
    Temp = 0;
    ExamBits (Temp, j, i, 0, comp.HighestBit, 1, Pow2, exam);
    return ;
  }

  Zone<MScal>::ZoneType No = Z->getNo ();
  MScal sup;
  sup = Z->getSup ();
  // if (q == 1)
  // q = 2;
  while (q <= sup) {
    if (timer.timeOver (config.duration))
      return ;
    if ((comp.implemCond == APP_FACT) && Z->DivQ[No])
      Quotient (comp.modulus.m, q, coef[j][i]);
    else
      coef[j][i] = q;
    ExamThisaj (j, i, Pow2, exam);
    if (Pow2)
      q += Eight;
    else
      ++q;
  }
}

//===========================================================================

// Builds random subregions of the research zones. This is only called for the
// Random search.
void ChoisirBornes (int j)
{
  long h;
  // Choisir une region au hasard et l'examiner au complet
  for (int i = 0; i < config.compon[j].k+1; i++) {
    if (i == config.compon[j].k)
      h = config.compon[j].Hk-1;
    else
      h = config.compon[j].H-1;
    Zone<MScal> *Z = zone[j] + i;
    Zone<MScal> *R = reg[j] + i;
    ++CoRegions[j];
    R->chooseBoundaries (config.compon[j], Z, h);
    // cout << (R + 1)->toString();
  }
}

//===========================================================================

// Exhaustive search in reg, a subset of zone, for the i-th coefficient of the
// j-th component. InsideExam will call this recursively if i < k.
void ExamRegion (int j, int i)
{
  // Used in the "Random search" case
  if (timer.timeOver (config.duration))
    return ;
  Zone<MScal> *R;
  R = i + reg[j];             // On va examiner toute cette region
  InsideExam (R, j, i, ExamRegion);
}

//===========================================================================

// Exhaustive search in the zone defined for the i-th coefficient of the j-th
// component. InsideExam will call this recursively if i < k.
void ExamAllZones (int j, int i)
{
  // This method is used in the {\em exhaustive search} case.
  if (timer.timeOver (config.duration))
    return ;
  Zone<MScal> *Z;
  Z = i + zone[j];

  while (Z != 0) {
    // On va examiner toute cette zone.
    InsideExam (Z, j, i, ExamAllZones);
    Z = Z->nextZone;
  }
}


//===========================================================================

// Allocates memory to local variables
void Init ()
{
  if (MINKL1 == config.normaType)
    Norm = L1NORM;

  pool = new LatticeTest<MScal, NScal> * [config.numGen[0]];
  // memset (pool, 0, config.numGen[0] * sizeof (LatticeTest *));

  coef.SetDims (config.J, config.getMaxK ()+1);

  CoNoElemPrim = new long[config.J];
  CoElemPrim = new long[config.J];
  TotEP = new long[config.J];
  CoPerMax = new long[config.J];
  CoRegions = new long[config.J];
  Hk = new long[config.J];
  H = new long[config.J];
  minVal = new double[1 + config.getMaxDim ()];
  for (int i = 0; i <= config.getMaxDim (); i++)
    minVal[i] = 0.0;

  compJ = new MRGComponent<MScal> * [config.J];
  primJ = new IntPrimitivity<MScal> * [config.J];

  int s;
  for (s = 0; s < config.J; s++) {
    Component<MScal> & comps = config.compon[s];
    CoNoElemPrim[s] = 0;
    CoElemPrim[s] = 0;
    TotEP[s] = 0;
    CoPerMax[s] = 0;
    CoRegions[s] = 0;

    cout << "1\n";
    if (comps.PerMax) {
      compJ[s] =
        new MRGComponent<MScal> (comps.modulus, comps.k,
            comps.F1, comps.file1.c_str (),
            comps.F2, comps.file2.c_str ());
      primJ[s] =
        new IntPrimitivity<MScal> (compJ[s]->ifm1, comps.modulus.m);
      cout << "2\n";

    } else
      compJ[s] = new MRGComponent<MScal> (comps.modulus.m, coef[s],  comps.k);
  }

  zone = new Zone<MScal> * [config.J];
  for (s = 0; s < config.J; s++)
    zone[s] = new Zone<MScal>[config.compon[s].k+1];

  reg = new Zone<MScal> * [config.J];
  for (s = 0; s < config.J; s++)
    reg[s] = new Zone<MScal>[config.compon[s].k+1];

  for (s = 0; s < config.J; s++) {
    for (int i = 0; i < config.compon[s].k; i++)
      coef[s][i] = 0;
  }
}


//===========================================================================

// Frees memory allocated to execute SeekMain
void Finalize ()
{
  int s;
  for (s = 0; s < config.J; s++)
    delete[]reg[s];
  delete[]reg;

  for (s = 0; s < config.J; s++)
    delete[]zone[s];
  delete[]zone;

  for (s = 0; s < config.J; s++) {
    if (config.compon[s].PerMax)
      delete primJ[s];
    delete compJ[s];
  }

  delete[]CoNoElemPrim;
  delete[]CoElemPrim;
  delete[]TotEP;
  delete[]CoPerMax;
  delete[]CoRegions;
  delete[]Hk;
  delete[]H;
  delete[]minVal;

  delete[]compJ;
  delete[]primJ;
  delete[]pool;
  delete normal;
}


//===========================================================================

// Initializes research zones. This fills the fields of the Zone objects in
// the zone array so that zones can then be used without errors.
void InitZones ()
{
  for (int s = 0; s < config.J; s++) {
    for (int i = 1; i < config.compon[s].k+1; i++) {
      zone[s][i].init (config.compon[s], s, i-1);
    }
  }
}


//===========================================================================

/* Examines the component j of the generator. Components from 0 to j-1 are 
 * already fixed. Subroutines called by this function will recursivelly call
 * SeekGen until all the components have been fixed. When this happens, SeekGen
 * is called with j = J and the current generator is tested
 * */
void SeekGen (int j)
{
  if (j >= config.J) {
    Test ();
      CompareMerit ();
      ++numATried;
      if(!(numATried%1000)) {
        std::cout << numATried << " generators tested.\n";
      }
      return ;
    }

  for (int region = 0; region < config.compon[j].numReg; region++) {
    if (timer.timeOver (config.duration))
      return ;
    if (config.compon[j].searchMethod == EXHAUST) {
      ExamAllZones (j, config.compon[j].k);
    } else {
      ChoisirBornes (j);
      ExamRegion (j, config.compon[j].k);
    }
  }
}

  //===========================================================================

  // Tests the generators that are in config.fileName. This simply performs the
  // tests as usual but does not go through the trouble of configuring
  // generators in a complicated way to work.
  void TestGen ()
  {
    int numGen, ln = 0;
    ParamReaderExt<MScal, NScal> reader(config.fileName);
    reader.getLines();
    reader.readInt(numGen, ln, 0);
    int k = config.compon[0].k;
    for (int i = 0; i < numGen; i++) {
      for (int j = 0; j < k; j++) {
        reader.readMScal(coef[0][j+1], ++ln, j);
      }
      Test ();
      CompareMerit ();
      ++numATried;
      if(!(numATried%1000)) {
        std::cout << numATried << " generators tested.\n";
      }
    }
    return ;
  }

}   // namespace


//===========================================================================

int main (int argc, char **argv)
{
  if (readConfigFile (argc, argv))
    return -1;
  //config.write ();
  Init ();
  timer.init ();
  Zone<MScal>::initFrontieres (config);
  InitZones ();
  if (config.readGenFile) {
    TestGen();
  } else {
    SeekGen (0);
  }
  SortBestGen ();
  if (config.outputType == RES) {
    PrintResults ();
  } else if (config.outputType == GEN) {
    PrintGen();
  } else {
    std::cout << "\nCould not print results as requested, tried to output using\
      standard method.\n";
    PrintResults ();
  }
  Finalize ();
}

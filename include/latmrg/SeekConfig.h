#ifndef SEEK_CONFIG_H
#define SEEK_CONFIG_H

#include <string>
#include <vector>

#include "latticetester/Types.h"
#include "latticetester/Const.h"

#include "latmrg/Modulus.h"
#include "latmrg/Const.h"


namespace LatMRG {

  /** This structure contains the information necessary to search an MRG 
   * component to a multiple MRG.
   * */
  template<typename Int>
    struct Component {
      /// Type of the generator. Either MRG or MWC.
      GenType genType;
      /// The modulus m of this component.
      Modulus<Int> modulus;
      /// The order of the component.
      int k;
      /** `bool` describing if the component has maximal period. It has to be 
       * `true` if the compenent has maximal period, and `false` otherwise.
       * */
      bool PerMax;
      /// \todo this
      ImplemCond implemCond;
      /// \todo this
      int NumBits;
      /// \todo this
      int HighestBit;
      /// \todo this
      int ncoef;
      /// \todo this
      int *Icoef;
      /// \todo this
      DecompType F1;
      /// \todo this
      std::string file1;
      /// \todo this
      DecompType F2;
      /// \todo this
      std::string file2;
      /// \todo this
      SearchMethod searchMethod;
      /// \todo this
      int numReg, H, Hk;
      /// \todo this
      MVect b, c;
      /// \todo this
      bool ApproxTotGen;
    };

  /**
   * \todo Compléter la description de la classe
   *
   * <strong> À EXAMINER:</strong> il pourrrait y avoir des avantages à inclure
   * une variable `LatConfig` à l’intérieur de `SeekConfig`. Plusieurs
   * variables n’auraient pas à être répétées dans `SeekConfig`, et la
   * recherche appelle les tests avec des paramètres de `LatConfig`, parfois
   * explicitement, ce qui cause des problèmes (cas PALPHA). Cela pourrait
   * simplifier le code.
   *
   */
  template<typename Int>
    class SeekConfig {
      public:
        static const int MAX_WORD_SIZE = 64;
        SeekConfig();
        ~SeekConfig();
        void write();
        Int* getMs() { return Ms; }
        int* getKs()  { return Ks; }
        int getMaxK();
        int getMaxTd();
        bool readGenFile; // read generator from file
        std::string fileName;  // file name
        int J;            // nombre de generateurs dans la combinason
        std::vector<Component<Int>> compon;
        Int* Ms;
        MVect* As;
        int* Ks;
        int C;
        double* minMerit;
        double* maxMerit; // C values
        int* numGen; // C values
        LatticeTester::CriterionType criter;
        LatticeTester::NormaType normaType;

        /**
         * The number of series of projections (see `td` below). The classical case
         * corresponds to \f$d=1\f$, for which the chosen test is done on all
         * successive dimensions from `td[0] = fromDim`, up to and including `td[1] =
         * toDim`.
         */
        // Criterion <Norm>

        int d;

        /**
         * Array containing the maximal dimensions for the projections in each
         * category. `td[1]` is the maximal dimension for successive
         * 1-dimensional projections, `td[2]` is the maximal dimension for
         * 2-dimensional projections, `td[3]` is the maximal dimension for
         * 3-dimensional projections, and so on. However, the value of `td[0]`
         * is the minimal dimension for the case of successive dimensions.
         */
        int *td;
        int getMaxDim() { return td[1]; }
        LatticeType latType;

        /**
         * If this flag is `true`, the test is to be applied on the dual lattice. If
         * it is `false`, the test is to be applied on the primal lattice.
         */
        bool dualF;

        /**
         * This flag applies specificaly to the spectral test by specifying at
         * what speed (and precision) it should be done. If it is set to `0` the
         * test is done normally without approximations. If it is set to `1`, the
         * test is done without Branch-and-Bound using only BKZ pre-reduction. 
         * Finaly, if it is set to `2` the test is done without Branch-and-Bound
         * using only LLL pre-reduction.
         * */
        int speed;

        /**
         * If `invertF` is `true`, the inverse of the length of the shortest
         * vector will be printed in the results. Otherwise, the length itself
         * will be printed.
         */
        bool invertF;

        /**
         * The value of \f$\alpha\f$ for the \f$P_{\alpha}\f$ test.
         */
        int alpha;
        int lacGroupSize;
        int lacSpacing;
        long maxNodesBB;
        double duration;
        long seed;  // seed of the random number generator
        LatticeTester::OutputType outputType;
    };

  //============================================================================

  template<typename Int>
    SeekConfig<Int>::SeekConfig ()
    {
      fileName.reserve(MAX_WORD_SIZE);
      td = new int[1 + d]; //
      C = 1;   // Number of categories; is 1 for now.
      minMerit = new double[C];
      maxMerit = new double[C];
      numGen = new int[C];
    }

  template<typename Int>
    SeekConfig<Int>::~SeekConfig()
    {
      fileName.clear();
      for (int i = 0;i < J;i++) {
        compon[i].file1.clear();
        compon[i].file2.clear();
        compon[i].b.kill();
        compon[i].c.kill();
      }
      delete[] td;
      delete[] minMerit;
      delete[] maxMerit;
      delete[] numGen;
    }

  template<typename Int>
    int SeekConfig<Int>::getMaxK()
    {
      int max = Ks[0];
      for (int i = 1; i < J; i++) {
        if (Ks[i] > max)
          max = Ks[i];
      }
      return max;
    }


  template<typename Int>
    int SeekConfig<Int>::getMaxTd()
    {
      int maxVal = 0;
      for (int i = 0; i < d; i++) {
        maxVal = std::max(maxVal, td[i]);
      }

      return maxVal;

    }

  template<typename Int>
    void SeekConfig<Int>::write()
    {
      std::cout << "readGenFile: " << std::boolalpha << readGenFile << std::endl;
      if (readGenFile)
        std::cout << "fileName: " << fileName << std::endl;
      std::cout << "J: " << J << std::endl;
      for (int i = 0; i < J; i++) {
        std::cout << "\n================ Component " << i+1 << " =================\n";
        std::cout << "typeGen: " << toStringGen (compon[i].genType) << std::endl;
        std::cout << "m: " << compon[i].modulus.m << std::endl;
        std::cout << "k: " << compon[i].k << std::endl;
        std::cout << "PerMax: " << compon[i].PerMax << std::endl;
        std::cout << "implemCond: " << toStringImplemCond (compon[i].implemCond) << std::endl;
        if (POWER_TWO == compon[i].implemCond) {
          std::cout << "NumBits: " << compon[i].NumBits << std::endl;
          std::cout << "HighestBit: " << compon[i].HighestBit << std::endl;
        }
        std::cout << "F1: " << toStringDecomp (compon[i].F1) << std::endl;
        std::cout << "file1: " << compon[i].file1 << std::endl;
        std::cout << "F2: " << toStringDecomp (compon[i].F2) << std::endl;
        std::cout << "file2: " << compon[i].file2 << std::endl;
        std::cout << "searchMethod: " << toStringSearchMethod (compon[i].searchMethod)
          << std::endl;
        if (RANDOM == compon[i].searchMethod) {
          std::cout << "numReg: " << compon[i].numReg << std::endl;
          std::cout << "H: " << compon[i].H << std::endl;
          std::cout << "Hk: " << compon[i].Hk << std::endl;
        }
        std::cout << "b: " << LatticeTester::toString<MVect>(compon[i].b, compon[i].k) << std::endl;
        std::cout << "c: " << LatticeTester::toString<MVect>(compon[i].c, compon[i].k) << std::endl;
      }
      std::cout << "========================= End  ======================\n" << std::endl;
      std::cout << "C : " << C << std::endl;
      std::cout << "MinMerit: " << LatticeTester::toString<double*>(minMerit, 0, C-1) << std::endl;
      std::cout << "MaxMerit: " << LatticeTester::toString<double*>(maxMerit, 0, C-1) << std::endl;
      std::cout << "NumGen: " << LatticeTester::toString<int*>(numGen, 0, C-1) << std::endl;
      std::cout << "td:   " << LatticeTester::toString (td, 0, d) << std::endl;
      std::cout << "criter: " << toStringCriterion(criter) << std::endl;
      std::cout << "normaType: " << toStringNorma (normaType) << std::endl;
      if (dualF)
        std::cout << "lattice: DUAL" << std::endl;
      else
        std::cout << "lattice: PRIMAL" << std::endl;
      std::cout << "latticeType: " << toStringLattice (latType) << std::endl;
      std::cout << "lacGroupSize: " << lacGroupSize << std::endl;
      std::cout << "lacSpacing: " << lacSpacing << std::endl;
      std::cout << "maxNodesBB: " << maxNodesBB << std::endl;
      std::cout << "duration: " << duration << std::endl;
      std::cout << "seed: " << seed << std::endl;
      //   std::cout << "S2: " << s2 << std::endl;
      std::cout << "outputType: " << LatticeTester::toStringOutput (outputType) << std::endl;
    }

}
#endif // SEEK_CONFIG_H

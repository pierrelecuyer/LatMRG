#ifndef LATMRG_LATTEST_H
#define LATMRG_LATTEST_H

template<typename Int, typename Dbl> struct LatTest {
  typedef NTL::vector<Int> IntVec;

  ConfigLat<Int, Dbl> conf;

  Chrono timer;
  int detail = 0;

  /**
   * Prints the results of the program execution.
   * */
  template<typename Lat>
    void printResults(MeritList<Lat>& bestLattice) {
      std::cout << "LatRe: A program to test Random Number Generators\n";
      std::cout << delim;
      std::cout << "Bellow are the results of the tests of " << conf.max_gen << " generators:\n";
      std::cout << "Generator type: " << toStringGen(conf.type[0]) << "\n";
      if (conf.type[0] == MRG) {
        std::cout << "With modulus:   m = " << conf.modulo[0] << " = " << conf.basis[0] << "^"
          << conf.exponent[0] << (conf.rest[0]>0?"+":"") << (conf.rest[0]!=0?std::to_string(conf.rest[0]):"") << "\n";
        std::cout << "Of order:       k = " << conf.order[0] << "\n";
      } else if (conf.type[0] == MWC) {
      } else if (conf.type[0] == MMRG) {
      }
      std::cout << "And " << (conf.period[0]?"full":"any") << " period length\n";
      std::cout << "The test was:\n";
      std::cout << "Minimal " << (conf.normaType==LatticeTester::NONE?"":"normalized ")
        << "shortest" << " non-zero vector length\n";
      if (conf.normaType != LatticeTester::NONE) {
        std::cout << "Normalizer used: "
          << LatticeTester::toStringNorma(conf.normaType) << "\n";
      }
      std::cout << "On dimensions and projections:\n";
      std::cout << conf.proj->toString();
      std::cout << delim;
      std::cout << "Allowed running time: " << conf.timeLimit << "s.\n";
      std::cout << "Actual CPU time: " << timer.toString() << "\n";
      std::cout << "Number of generators tested: " << conf.max_gen << "\n\n";
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        std::cout << delim;
        if (conf.type[0] == MRG) {
          std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
        } else if (conf.type[0] == MWC) {}
        else if (conf.type[0] == MMRG) {}
        else if (conf.type[0] == COMBO) {
          std::cout << (*it).getLattice();
        }
        if (detail == 0) {
          std::cout << (*it).toStringMerit();
        } else if (detail == 1) {
          std::cout << (*it).toStringDim();
        } else if (detail == 2) {
          std::cout << (*it).toStringProjections();
        }
      }
    }

  //==========================================================================

  int TestLat () {
    if (!conf.gen_set) {
      std::cerr << "No generator set for in seek tag. Aborting.\n";
      return 1;
    }
    if (!(conf.test_set)) {
      std::cerr << "No test set for in seek tag. Aborting.\n";
      return 1;
    }
    if (!conf.proj_set) {
      std::cerr << "No projections set for in seek tag. Aborting.\n";
      return 1;
    }
    // Initializing values
    timer.init();

    // Testing the generator(s)
    if (conf.num_comp > 1) {
      // Combined generators case
      MeritList<ComboLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MRGLattice<Int, Dbl>* mrg = getLatCombo<Int, Dbl>(conf.fact, conf.max_dim);
      ComboLattice<Int, Dbl> combolat(conf.fact, *mrg);
      bestLattice.add(test_lat(combolat, conf));
      printResults(bestLattice);
    } else if (conf.type[0] == MRG) {
      MeritList<MRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      IntVec temp(conf.order[0]+1);
      temp[0] = Int(0);
      for (int i = 1; i < conf.order[0]+1; i++) temp[i] = conf.coeff[0][i-1];
      MRGLattice<Int, Dbl> mrglat(conf.modulo[0], temp, conf.max_dim, conf.order[0], FULL);
      bestLattice.add(test_lat(mrglat, conf));
      printResults(bestLattice);
    } else if (conf.type[0] == MWC) {
      //MWCLattice<Int, Dbl>* mwclat = 0;
      //mwclat = nextGenerator(mwclat);
      //Dbl merit(test(*mwclat));
      //if (merit > bestMerit) {
      //  addLattice(mwclat, merit);
      //}
    } else if (conf.type[0] == MMRG) {
      MeritList<MMRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MMRGLattice<Int, Dbl> mmrglat(conf.modulo[0], conf.matrix, conf.max_dim, conf.order[0]);
      bestLattice.add(test_lat(mmrglat, conf));
      printResults(bestLattice);
    }
    delete conf.proj;
    return 0;
  }

}; // end struct LatTest
#endif

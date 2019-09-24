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
      std::cout << "Generator type: " << toStringGen(conf.fact[0]->get_type()) << "\n";
      if (conf.fact[0]->get_type() == MRG) {
        std::cout << "With modulus:   m = " << conf.fact[0]->getM() << " = " << conf.fact[0]->getB() << "^"
          << conf.fact[0]->getE();
       if (conf.fact[0]->getR() > 0) std::cout << "+" << conf.fact[0]->getR();
       std::cout << "\n";
        std::cout << "Of order:       k = " << conf.fact[0]->getK() << "\n";
      } else if (conf.fact[0]->get_type() == MWC) {
      } else if (conf.fact[0]->get_type() == MMRG) {
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
        if (conf.fact[0]->get_type() == MRG) {
          std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
        } else if (conf.fact[0]->get_type() == MWC) {}
        else if (conf.fact[0]->get_type() == MMRG) {}
        else if (conf.fact[0]->get_type() == COMBO) {
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
    } else if (conf.fact[0]->get_type() == MRG) {
      MeritList<MRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      IntVec temp(conf.fact[0]->getK()+1);
      temp[0] = Int(0);
      for (int i = 1; i < conf.fact[0]->getK()+1; i++) temp[i] = conf.coeff[0][i-1];
      MRGLattice<Int, Dbl> mrglat(conf.fact[0]->getM(), temp, conf.max_dim, conf.fact[0]->getK(), FULL);
      bestLattice.add(test_lat(mrglat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MWC) {
      //MWCLattice<Int, Dbl>* mwclat = 0;
      //mwclat = nextGenerator(mwclat);
      //Dbl merit(test(*mwclat));
      //if (merit > bestMerit) {
      //  addLattice(mwclat, merit);
      //}
    } else if (conf.fact[0]->get_type() == MMRG) {
      MeritList<MMRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MMRGLattice<Int, Dbl> mmrglat(conf.fact[0]->getM(), conf.matrix, conf.max_dim, conf.fact[0]->getK());
      bestLattice.add(test_lat(mmrglat, conf));
      printResults(bestLattice);
    }
    delete conf.proj;
    for (int i = 0; i < conf.num_comp; i++) delete conf.fact[i];
    return 0;
  }

}; // end struct LatTest
#endif

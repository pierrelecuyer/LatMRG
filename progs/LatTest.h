#ifndef LATMRG_LATTEST_H
#define LATMRG_LATTEST_H

template<typename Int, typename Dbl> struct LatTest {
  typedef NTL::vector<Int> IntVec;

  ConfigLat<Int, Dbl> conf;

  Chrono timer;
  int detail = 0;
  std::vector<bool> full_period;

  /**
   * Prints the results of the program execution.
   * */
  template<typename Lat>
    void printResults(MeritList<Lat>& bestLattice) {
      std::cout << "LatTest: A program to test Random Number Generators\n";
      std::cout << delim;
      std::cout << ((conf.num_comp>1)?"Combined generators":"Simple generator")
        << " configuration" << ((conf.num_comp>1)?"s":"") << "\n\n";
      for (int k = 0; k < conf.num_comp; k++) {
        if (k > 0) std::cout << "\n";
        if (conf.num_comp >1) std::cout << "Component " << k+1 << ":\n";
        std::cout << "Generator type: " << toStringGen(conf.fact[k]->get_type()) << "\n";
        if (conf.fact[k]->get_type() == MRG) {
          std::cout << "Modulus:        m = " << conf.fact[k]->getM() << " = " << conf.fact[k]->getB() << "^"
            << conf.fact[k]->getE();
          if (conf.fact[k]->getR() > 0) std::cout << "+" << conf.fact[k]->getR();
          if (conf.fact[k]->getR() < 0) std::cout << conf.fact[k]->getR();
          std::cout << "\n";
          std::cout << "Order:          k = " << conf.fact[k]->getK() << "\n";
        } else if (conf.fact[k]->get_type() == MWC) {
        } else if (conf.fact[k]->get_type() == MMRG) {
        }
        std::cout << (conf.period[0]?"Check":"Don't check") << " full period length\n";
      }
      std::cout << "\nTest:\n";
      if (conf.criterion == LatticeTester::SPECTRAL) {
        std::cout << "Spectral Test\n";
        if (conf.normaType != LatticeTester::NONE) std::cout << "Normalizer used: "
          << toStringNorma(conf.normaType);
      } else if (conf.criterion == LatticeTester::BEYER) std::cout << "Beyer quotient";
      else if (conf.criterion == LatticeTester::LENGTH) std::cout << "Shortest vector length";
      std::cout << "\n\n";
      std::cout << "Dimensions and projections:\n";
      std::cout << conf.proj->toString();
      std::cout << delim;
      std::cout << "Allowed running time: " << conf.timeLimit << "s.\n";
      std::cout << "Actual CPU time: " << timer.toString() << "\n";
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        std::cout << delim;
        std::cout << (*it).getLattice();
        if (conf.num_comp > 1) {
          for (int i = 0; i<conf.num_comp; i++) {
            if (conf.period[i]) std::cout << "Component " << i+1
              << ((full_period[i])?" has":" does not have") << " full period.\n";
          }
        } else {
          if (conf.period[0]) std::cout << "Full period: " << full_period[0] << "\n";
        }
        std::cout << "\n";
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
      full_period.resize(conf.num_comp);
      // Checking full period of components that require it
      for (int i = 0; i < conf.num_comp; i++) {
        IntVec temp(conf.fact[i]->getK()+1);
        temp[0] = Int(0);
        for (int j = 1; j < conf.fact[i]->getK()+1; j++) temp[j] = conf.coeff[i][j-1];
        if (conf.period[i]) full_period[i] = conf.fact[i]->maxPeriod(temp);
        // This is not set because it is not needed in Seek
        conf.fact[i]->setA(conf.coeff[i]);
      }
      // Combined generators case
      MeritList<ComboLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MRGLattice<Int, Dbl>* mrg = getLatCombo<Int, Dbl>(conf.fact, conf.max_dim);
      ComboLattice<Int, Dbl> combolat(conf.fact, *mrg);
      bestLattice.add(test_lat(combolat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MRG) {
      full_period.resize(1);
      MeritList<MRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      IntVec temp(conf.fact[0]->getK()+1);
      temp[0] = Int(0);
      for (int i = 1; i < conf.fact[0]->getK()+1; i++) temp[i] = conf.coeff[0][i-1];
      if (conf.period[0]) full_period[0] = conf.fact[0]->maxPeriod(temp);
      MRGLattice<Int, Dbl> mrglat(conf.fact[0]->getM(), temp, conf.max_dim, conf.fact[0]->getK(), FULL);
      bestLattice.add(test_lat(mrglat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MWC) {
      full_period.resize(1);
      //MWCLattice<Int, Dbl>* mwclat = 0;
      //mwclat = nextGenerator(mwclat);
      //Dbl merit(test(*mwclat));
      //if (merit > bestMerit) {
      //  addLattice(mwclat, merit);
      //}
    } else if (conf.fact[0]->get_type() == MMRG) {
      full_period.resize(1);
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

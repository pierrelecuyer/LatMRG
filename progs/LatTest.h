#ifndef LATMRG_LATTEST_H
#define LATMRG_LATTEST_H

extern std::ostream* out;

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
      *out << "LatTest: A program to test Random Number Generators\n";
      *out << delim;
      *out << ((conf.num_comp>1)?"Combined generators":"Simple generator")
        << " configuration" << ((conf.num_comp>1)?"s":"") << "\n\n";
      for (int k = 0; k < conf.num_comp; k++) {
        if (k > 0) *out << "\n";
        if (conf.num_comp >1) *out << "Component " << k+1 << ":\n";
        *out << "Generator type: " << toStringGen(conf.fact[k]->get_type()) << "\n";
        *out << "Modulo:         m = " ;
        if (conf.fact[k]->get_type() == MWC) *out << conf.fact[k]->m_MWCb;
        else *out << conf.fact[k]->getM();
        *out << " = " << conf.fact[k]->getB() << "^" << conf.fact[k]->getE();
        if (conf.fact[k]->getR() > 0) *out << "+" << conf.fact[k]->getR();
        if (conf.fact[k]->getR() < 0) *out << conf.fact[k]->getR();
        *out << "\n";
        *out << "Order:          k = " << conf.fact[k]->getK() << "\n";
        *out << (conf.period[0]?"Check":"Don't check") << " full period length\n";
      }
      *out << "\nTest:\n";
      if (conf.criterion == LatticeTester::SPECTRAL) {
        *out << "Spectral Test";
        if (conf.normaType != LatticeTester::NONE) *out << "\nNormalizer used: "
          << toStringNorma(conf.normaType);
      } else if (conf.criterion == LatticeTester::BEYER) *out << "Beyer quotient";
      else if (conf.criterion == LatticeTester::LENGTH) *out << "Shortest vector length";
      *out << "\n\n";
      *out << "Dimensions and projections:\n";
      *out << conf.proj->toString();
      *out << delim;
      *out << "Allowed running time: " << conf.timeLimit << "s.\n";
      *out << "Actual CPU time: " << timer.toString() << "\n";
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        *out << delim;
        *out << (*it).getLattice() << "\n";
        if (conf.num_comp > 1) {
          bool print = false;
          for (int i = 0; i<conf.num_comp; i++) {
            if (conf.period[i]){
              *out << "Component " << i+1
              << ((full_period[i])?" has":" does not have") << " full period.\n";
              print = true;
            }
          }
          if (print) *out << "\n";
        } else {
          if (conf.period[0]) *out << "Full period: " << (full_period[0]?"yes":"no") << "\n\n";
        }
        if (detail == 0) {
          *out << (*it).toStringMerit();
        } else if (detail == 1) {
          *out << (*it).toStringDim();
        } else if (detail == 2) {
          *out << (*it).toStringProjections();
        }
      }
    }

  //==========================================================================

  int TestLat () {
    if (!conf.gen_set) {
      std::cerr << "No generator set for in lattest tag. Aborting.\n";
      return 1;
    }
    if (!(conf.test_set)) {
      std::cerr << "No test set for in lattest tag. Aborting.\n";
      return 1;
    }
    if (!conf.proj_set) {
      std::cerr << "No projections set for in lattest tag. Aborting.\n";
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
      MeritList<MWCLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MWCLattice<Int, Dbl> mwclat(conf.fact[0]->m_MWCb, conf.fact[0]->getM());
      if (conf.period[0]) full_period[0] = conf.fact[0]->maxPeriod(mwclat.getCoef());
      bestLattice.add(test_lat(mwclat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MMRG) {
      full_period.resize(1);
      if (conf.period[0]) full_period[0] = conf.fact[0]->maxPeriod(conf.fact[0]->getMatrix());
      MeritList<MMRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      MMRGLattice<Int, Dbl> mmrglat(conf.fact[0]->getM(), conf.fact[0]->getMatrix(), conf.max_dim, conf.fact[0]->getK());
      bestLattice.add(test_lat(mmrglat, conf));
      printResults(bestLattice);
    }
    delete conf.proj;
    for (int i = 0; i < conf.num_comp; i++) delete conf.fact[i];
    return 0;
  }

}; // end struct LatTest
#endif

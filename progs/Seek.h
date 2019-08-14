#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

// Number of generators to generate before testing time in nextGenerator
#define DELAY 1000

template<typename Int, typename Dbl> struct SeekMain {
  typedef NTL::vector<Int> IntVec;
  typedef NTL::matrix<Int> IntMat;

  ConfigSeek<Int, Dbl> conf;

  // Program global objects
  Chrono timer; // program timer
  int maxDim = 40;

  // MRG specific parameters
  MRGComponent<Int>* mrg;
  int** comb_coeff;

  /**
   * A small class to search for modulus for MWC generators.
   * */
  class Modulus {
    public:

      /**
       * Search will begin at `2^e±c` and will increase of decrease depending on
       * `increase`.
       * */
      Modulus(long e, Int c, bool increase) {
        m = (Int(1)<<e) - 1;
        this->increase = increase;
        if (increase) m+=c;
        else m-=c;
      }

      /**
       * This function finds the next value for `m`. This can return 1 if
       * `increase` is `false` and `m` gets to 1.
       * */
      Int next() {
        while (m > 1) {
          nextM();
          LatticeTester::PrimeType status = LatticeTester::IntFactor<Int>::isPrime (m, KTRIALS);
          if (status == LatticeTester::PRIME || status == LatticeTester::PROB_PRIME) {
            if (1 == m % 4) continue;
            Int m1 = (m - 1)/2;
            status = LatticeTester::IntFactor<Int>::isPrime (m1, KTRIALS);
            if (status != LatticeTester::PRIME && status != LatticeTester::PROB_PRIME) continue;
            // For MWC we check 2^{(m-1)/2} \neq 1 mod m.
            NTL::ZZ_p::init(NTL::ZZ(m));
            NTL::ZZ_p a = NTL::power(NTL::ZZ_p(2), m1);
            if (NTL::IsOne(a)) {
              continue;
            }
            return m;
          }
        }
        return Int(-1);
      }

    private:

      static const long KTRIALS = 200;

      /**
       * Last valid modulus found (or 2^e).
       * */
      Int m;

      /**
       * As passed to constructor.
       * */
      bool increase;

      /**
       * Increment/decrement m
       * */
      void nextM () {
        if (increase) m += 2;
        else m -= 2;
        if (0 == m % 5){
          if (increase) m += 2;
          else m -= 2;
        }
      }
  } /**mod*/;

  /*
   * The goal is to create this overload and to use it to switch generators
   * without requiring the use of switch statements.
   * */
  MRGLattice<Int, Dbl>* nextGenerator(MRGLattice<Int, Dbl>* lattice) {
    // Setting up two vectors. MRGComponent and MRGLattice do not use the same
    // vector format
    IntVec A;
    A.SetLength(conf.order+1);
    NTL::clear(A);
    int delay = 0;
    // The program will not run the maxPeriod function if it is not wanted with
    // this condition
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(conf.timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<conf.order; i++) A[i+1] = conf.coeff[i] * randInt(Int(0), conf.modulo);
      delay++;
    } while ((A[conf.order] == 0) || (conf.period && !mrg->maxPeriod(A)));
    if (lattice) delete lattice;
    return new MRGLattice<Int, Dbl>(conf.modulo, A, maxDim, conf.order, FULL);
  }

  MRGLattice<Int, Dbl>* nextGeneratorPow2(MRGLattice<Int, Dbl>* lattice) {
    // Setting up two vectors. MRGComponent and MRGLattice do not use the same
    // vector format
    IntVec A;
    A.SetLength(conf.order+1);
    NTL::clear(A);
    IntVec coefficients(2*conf.order);
    int sign;
    int delay = 0;
    // The program will not run the maxPeriod function if it is not wanted with
    // this condition
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(conf.timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<conf.order; i++) {
        if (conf.coeff[2*i] < 0) {
          // This is a placeholder value for a zero coefficient
          coefficients[2*i] = coefficients[2*i+1] = 2004012;
          A[i+1] = 0;
          continue;
        }
        coefficients[2*i] = randInt(Int(LatticeTester::Lg(conf.rest))+1, conf.coeff[2*i]);
        sign = randInt(0,1);
        {
          Int tmp;
          NTL::power2(tmp, coefficients[2*i]);
          A[i+1] = Int(sign?1:-1) * tmp;
        }
        coefficients[2*i] ^= sign<<30;
        if (!(conf.coeff[2*i+1] < 0)) {
          coefficients[2*i+1] = randInt(Int(LatticeTester::Lg(conf.rest))+1, conf.coeff[2*i+1]);
          sign = randInt(0,1);
          Int tmp;
          NTL::power2(tmp, coefficients[2*i+1]);
          A[i+1] += Int(sign?1:-1) * tmp;
          coefficients[2*i+1] ^= (sign<<30);
        }
        else coefficients[2*i+1] = 2004012;
      }
      delay++;
    } while ((A[conf.order] == 0) || (conf.period && !mrg->maxPeriod(A)));
    if (lattice) delete lattice;
    MRGLattice<Int, Dbl>* lat = new MRGLattice<Int, Dbl>(conf.modulo, A, maxDim, conf.order, FULL);
    lat->setPower2(coefficients);
    return lat;
  }

  MWCLattice<Int, Dbl>* nextGenerator(MWCLattice<Int, Dbl>* lattice) {
    Int m(0);
    long exp = conf.exponent-1;
    // 63 bits at a time because NTL converts from SIGNED long
    while(exp > 0) {
      if (exp < 63) {
        m = (m << exp) + LatticeTester::RandBits(exp);
        exp -= exp;
      }
      m = (m << 63) + LatticeTester::RandBits(63);
      exp -= 63;
    }
    if ((m&1) == 1) m+=1;
    Modulus mod(conf.exponent, m, true);
    if (lattice) delete lattice;
    return new MWCLattice<Int, Dbl>(conf.b, mod.next());
  }

  /*
   * Goin' full random for now
   * */
  MMRGLattice<Int, Dbl>* nextGenerator(MMRGLattice<Int, Dbl>* lattice) {
    IntMat A;
    A.SetDims(conf.order, conf.order);
    NTL::clear(A);
    int delay = 0;
    do {
      if (delay >= DELAY) {
        if (timer.timeOver(conf.timeLimit)) return NULL;
        else delay = 0;
      }
      for (long i = 0; i<conf.order-1; i++) {
        A[i][i+1] = Int(1);
      }
      for (long i = 0; i<conf.order; i++) {
        A[conf.order-1][i] = randInt(Int(0), conf.modulo);
      }
      delay++;
    } while ((NTL::determinant(A) == 0) || (conf.period && !mrg->maxPeriod(A)));
    // Correcting the matrix to a full matrix
    //for (int i = 0; i<conf.order; i++) A *= A;
    for (int i = 0; i<conf.order; i++)
      for (int j = 0; j<conf.order; j++)
        A[i][j] = A[i][j]%conf.modulo;
    if (lattice) delete lattice;
    return new MMRGLattice<Int, Dbl>(conf.modulo, A, maxDim, conf.order);
  }

  ComboLattice<Int, Dbl>* nextGenerator(ComboLattice<Int, Dbl>* lattice) {
    // Setting up two vectors. MRGComponent and MRGLattice do not use the same
    // vector format
    std::vector<MRGComponent<Int>> components;
    for (int k = 0; k < conf.num_comp; k++) {
      IntVec A;
      A.SetLength(conf.comb_order[k]+1);
      NTL::clear(A);
      IntVec coefficients(2*conf.comb_order[k]);
      int sign;
      int delay = 0;
      // The program will not run the maxPeriod function if it is not wanted with
      // this condition
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.comb_order[k]; i++) {
          if (comb_coeff[k][2*i] < 0) {
            // This is a placeholder value for a zero coefficient
            coefficients[2*i] = coefficients[2*i+1] = 2004012;
            A[i+1] = 0;
            continue;
          }
          coefficients[2*i] = randInt(int(LatticeTester::Lg(abs(conf.comb_rest[k])))+1, comb_coeff[k][2*i]);
          sign = randInt(0,1);
          {
            Int tmp;
            NTL::power2(tmp, coefficients[2*i]);
            A[i+1] = Int(sign?1:-1) * tmp;
          }
          coefficients[2*i] ^= sign<<30;
          if (!(comb_coeff[k][2*i+1] < 0)) {
            coefficients[2*i+1] = randInt(int(LatticeTester::Lg(abs(conf.comb_rest[k])))+1, comb_coeff[k][2*i+1]);
            sign = randInt(0,1);
            Int tmp;
            NTL::power2(tmp, coefficients[2*i+1]);
            A[i+1] += Int(sign?1:-1) * tmp;
            coefficients[2*i+1] ^= (sign<<30);
          }
          else coefficients[2*i+1] = 2004012;
        }
        delay++;
      } while ((A[conf.comb_order[k]] == 0) || (conf.comb_period[k] && !conf.comb_fact[k]->maxPeriod(A)));
      IntVec B;
      B.SetLength(conf.comb_order[k]);
      for (int i = 0; i < conf.comb_order[k]; i++) B[i] = A[i+1];
      components.push_back(MRGComponent<Int>(conf.comb_modulo[k], B, conf.comb_order[k]));
    }
    if (lattice) delete lattice;
    MRGLattice<Int, Dbl>* mrg_lat = getLatCombo<Int, Dbl>(components, maxDim);
    ComboLattice<Int, Dbl>* new_lat = new ComboLattice<Int, Dbl>(components, *mrg_lat);
    delete mrg_lat;
    return new_lat;
  }

  /*
   * These next function add the tested lattices to the list of the best ones.
   * This only add the lattices that are good enough.
   * */
  template<typename Lat>
    void printResults(MeritList<Lat>& bestLattice) {
      std::cout << "\nSeekRe: A search program for Random Number Generators\n";
      std::cout << delim;
      std::cout << "Bellow are the results of a search for random number generators:\n";
      std::cout << "Generator type: " << toStringGen(conf.type) << "\n";
      if (conf.type == MRG) {
        std::cout << "With modulus:   m = " << conf.modulo << " = " << conf.basis << "^"
          << conf.exponent << (conf.rest>0?"+":"") << (conf.rest!=0?std::to_string(conf.rest):"") << "\n";
        std::cout << "Of order:       k = " << conf.order << "\n";
      } else if (conf.type == MWC) {
      } else if (conf.type == MMRG) {
      }
      std::cout << "And " << (conf.period?"full":"any") << " period length\n";
      std::cout << "The test was:\n" << (conf.best?"Best":"Worst") << " generators "
        "ranked by ";
      if(conf.criterion == LatticeTester::SPECTRAL) std::cout << "minimal "
        << (conf.normaType==LatticeTester::NONE?"inverse":"normalized")
          << " shortest non-zero vector length (Spectral test)\n";
      else if (conf.criterion == LatticeTester::LENGTH) std::cout << "minimal"
        << " shortest non-zero vector length\n";
      else if (conf.criterion == LatticeTester::BEYER) std::cout << "their Beyer quotient\n";
      if (conf.normaType != LatticeTester::NONE) {
        std::cout << "Normalizer used: "
          << LatticeTester::toStringNorma(conf.normaType) << "\n";
      }
      std::cout << "On dimensions and projections:\n";
      std::cout << conf.proj->toString();
      std::cout << delim;
      std::cout << "Allowed running time: " << conf.timeLimit << "s.\n";
      std::cout << "Actual CPU time: " << timer.toString() << "\n";
      std::cout << "Number of generators kept: " << conf.max_gen << "\n";
      std::cout << "Number of generators tested: " << conf.num_gen << "\n\n";
      std::cout << "Retained generators (from best to worst):\n";
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        std::cout << delim;
        if (conf.type == MRG) {
          std::cout << "Coefficients:\n" << (*it).getLattice() << "\n";
        } else if (conf.type == MWC) {}
        else if (conf.type == MMRG) {
          std::cout << "Matrix:\n" << (*it).getLattice() << "\n";
        } else if (conf.type == COMBO) {
          std::cout << (*it).getLattice() << "\n";
        }
        std::cout << (*it).toStringMerit();
      }
    }

  int print_progress(int old) {
    int per_80 = timer.val(Chrono::SEC)/conf.timeLimit * 80;
    if (per_80 > 80) per_80 = 80;
    if (per_80 < 0) per_80 = 0;
    // We do not print for no reason as this slows the program a lot.
    if (per_80 <= old) return old;
    std::cout << "[";
    for (int i = 0; i < per_80; i++) std::cout << "#";
    for (int i = per_80; i < 80; i++) std::cout << " ";
    std::cout << "] ";
    std::cout << std::setw(2) << int(per_80/80.0*100) << " %\r" << std::flush;
    return per_80;
  }

  int Seek ()
  {
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
    srand(time(NULL));
    conf.filem1 = "./tempm1" + std::to_string(rand());
    conf.filer = "./tempr" + std::to_string(rand());
    // Dynamically allocated objects
    if (conf.type != COMBO) {
      mrg = new MRGComponent<Int>(conf.modulo, conf.order, conf.decompm1,
          conf.filem1.c_str(), conf.decompr, conf.filer.c_str());
    }
    // std::vector<int> projDim = {31};
    // conf.proj = new Projections(1, 7, projDim);
    timer.init();

    // Launching the tests
    std::cout << "Program progress:\n";
    int old = print_progress(-1);
    if (conf.type == MRG) {
      MRGLattice<Int, Dbl>* mrglat = 0;
      MeritList<MRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      while (!timer.timeOver(conf.timeLimit)) {
        if (conf.construction == "POW2") mrglat = nextGeneratorPow2(mrglat);
        else if (conf.construction == "RANDOM") mrglat = nextGenerator(mrglat);
        if (mrglat == NULL) continue;
        bestLattice.add(test_seek(*mrglat, conf));
        conf.num_gen++;
        conf.currentMerit = bestLattice.getMerit();
        old = print_progress(old);
      }
      if (mrglat) delete mrglat;
      printResults(bestLattice);
    } else if (conf.type == MWC) {
      MWCLattice<Int, Dbl>* mwclat = 0;
      MeritList<MWCLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      while (!timer.timeOver(conf.timeLimit)) {
        mwclat = nextGenerator(mwclat);
        if (mwclat == NULL) continue;
        bestLattice.add(test_seek(*mwclat, conf));
        conf.num_gen++;
        conf.currentMerit = bestLattice.getMerit();
        old = print_progress(old);
      }
      if (mwclat) delete mwclat;
      printResults(bestLattice);
    } else if (conf.type == MMRG) {
      MMRGLattice<Int, Dbl>* mmrglat = 0;
      MeritList<MMRGLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      while (!timer.timeOver(conf.timeLimit)) {
        mmrglat = nextGenerator(mmrglat);
        if (mmrglat == NULL) continue;
        bestLattice.add(test_seek(*mmrglat, conf));
        conf.num_gen++;
        conf.currentMerit = bestLattice.getMerit();
        old = print_progress(old);
      }
      if (mmrglat) delete mmrglat;
      printResults(bestLattice);
    } else if (conf.type == COMBO) {
      ComboLattice<Int, Dbl>* combolat=0;
      MeritList<ComboLattice<Int, Dbl>> bestLattice(conf.max_gen, true);
      while (!timer.timeOver(conf.timeLimit)) {
        combolat = nextGenerator(combolat);
        if (combolat == NULL) continue;
        bestLattice.add(test_seek(*combolat, conf));
        conf.num_gen++;
        conf.currentMerit = bestLattice.getMerit();
        old = print_progress(old);
      }
      if (combolat) delete combolat;
      printResults(bestLattice);
      for (int i = 0; i < conf.num_comp; i++) {
        if (conf.comb_fact[i]) delete conf.comb_fact[i];
        delete[] comb_coeff[i];
      }
      delete[] comb_coeff;
    }
    delete conf.proj;
    delete mrg;
    return 0;
  }
}; // end struct SeekMain
#endif

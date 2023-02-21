#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

// Number of generators to generate before testing time in nextGenerator
#define DELAY 1000
extern std::ostream* out;
extern bool print_time;

namespace {
  // Program global objects
  Chrono timer; // program timer
}

/**
 * A namespace containing all the search functions.
 * */
namespace LatMRGSeek {
  /**
   * A small class to search for modulus for MWC generators.
   * */
  template<typename Int>
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
  // extern template class Modulus<std::int64_t>;
  // extern template class Modulus<NTL::ZZ>;

  /*
   * The goal is to create this overload and to use it to switch generators
   * without requiring the use of switch statements.
   * */
  template<typename Int, typename Real>
    MRGLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf) {
      // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
      // vector format
      NTL::vector<Int> A;
      A.SetLength(conf.fact[0]->getK()+1);
      NTL::clear(A);
      int delay = 0;
      // The program will not run the maxPeriod function if it is not wanted with
      // this condition
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.fact[0]->getK(); i++) A[i+1] = randInt(Int(0), conf.coeff[0][i]);
        delay++;
      } while ((A[conf.fact[0]->getK()] == 0) ||
          (conf.period[0] && (conf.fact[0]->get_type()==MRG) ?
           !conf.fact[0]->maxPeriod(A) :
           conf.fact[0]->maxPeriod(A[conf.fact[0]->getK()])));
      return new MRGLattice<Int, Real>(conf.fact[0]->getM(), A, conf.proj->numProj(), conf.fact[0]->getK(), FULL);
    }

  /*
   * The goal is to create this overload and to use it to switch generators
   * without requiring the use of switch statements.
   * */
  template<typename Int, typename Real>
    MRGLattice<Int, Real>* MRGApproxFactor(ConfigSeek<Int, Real>& conf) {
      // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
      // vector format
      NTL::vector<Int> A;
      A.SetLength(conf.fact[0]->getK()+1);
      NTL::clear(A);
      int delay = 0;
      // The program will not run the maxPeriod function if it is not wanted with
      // this condition
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.fact[0]->getK(); i++) A[i+1] = conf.coeff[0][i] * randInt(Int(0), NTL::SqrRoot(conf.fact[0]->getM()));
        delay++;
      } while ((A[conf.fact[0]->getK()] == 0) ||
          (conf.period[0] && (conf.fact[0]->get_type()==MRG) ?
           !conf.fact[0]->maxPeriod(A) :
           conf.fact[0]->maxPeriod(A[conf.fact[0]->getK()])));
      return new MRGLattice<Int, Real>(conf.fact[0]->getM(), A, conf.proj->numProj(), conf.fact[0]->getK(), FULL);
    }

  /**
   * Enuemrates all the possible generators for this generator. Needs to
   * remember previous generator and last modified coordinate between calls.
   * */
  template<typename Int, typename Real>
    class MRGLatticeExhaust {
      static NTL::vector<Int> A;
      public:
        static MRGLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf) {
          // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
          // vector format
          if (conf.num_gen == 0) {
            A.SetLength(conf.fact[0]->getK()+1);
          }
          //int delay = 0;
          // The program will not run the maxPeriod function if it is not wanted with
          // this condition
          do {
            // if (delay >= DELAY) {
            //   if (timer.timeOver(conf.timeLimit)) return NULL;
            //   else delay = 0;
            // }
            for (int i = 1; i <= conf.fact[0]->getK(); i++) {
              if (conf.coeff[0][i-1] == 0) continue;
              A[i]++;
              if (i == conf.fact[0]->getK()) {
              }
              if (A[i] < conf.fact[0]->getM()) break;
              else {
                if (i == conf.fact[0]->getK()) {
                  return NULL;
                }
                A[i] = 0;
              }
            }
            //delay++;
            // std::cout << "A " << A << "\n";
          } while ((A[conf.fact[0]->getK()] == 0) ||
              (conf.period[0] && (conf.fact[0]->get_type()==MRG) ?
               !conf.fact[0]->maxPeriod(A) :
               conf.fact[0]->maxPeriod(A[conf.fact[0]->getK()])));
          return new MRGLattice<Int, Real>(conf.fact[0]->getM(), A, conf.proj->numProj(), conf.fact[0]->getK(), FULL);
        }
    };
  template<> NTL::vector<std::int64_t> MRGLatticeExhaust<std::int64_t, double>::A(0);
  template<> NTL::vector<NTL::ZZ> MRGLatticeExhaust<NTL::ZZ, double>::A(0);
  template<> NTL::vector<NTL::ZZ> MRGLatticeExhaust<NTL::ZZ, NTL::RR>::A(0);

  template<typename Int, typename Real>
    MRGLattice<Int, Real>* nextGeneratorPow2(ConfigSeek<Int, Real>& conf) {
      typedef NTL::vector<Int> IntVec;
      // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
      // vector format
      IntVec A;
      A.SetLength(conf.fact[0]->getK()+1);
      int delay = 0;
      std::vector<IntVec> coeffs;
      coeffs.resize(conf.fact[0]->getK());
      for (int i = 0; i < conf.fact[0]->getK(); i++) coeffs[i].resize(
          NTL::conv<int>(conf.coeff[0][conf.fact[0]->getK()]));
      // The program will not run the maxPeriod function if it is not wanted with
      // this condition
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.fact[0]->getK(); i++) {
          A[i+1] = 0;
          if (conf.coeff[0][i] >= 0) {
            for (long j = 0; j < conf.coeff[0][conf.fact[0]->getK()]; j++) {
              coeffs[i][j] = randInt(Int(0), conf.coeff[0][i]);
              Int tmp2;
              NTL::power2(tmp2, coeffs[i][j]);
              A[i+1] += Int(randInt(0,1)?-1:1) * tmp2;
            }
          }
        }
        delay++;
      } while ((A[conf.fact[0]->getK()] == 0) || (conf.period[0] && !conf.fact[0]->maxPeriod(A)));
      MRGLattice<Int, Real>* lat = new MRGLattice<Int, Real>(conf.fact[0]->getM(), A, conf.proj->numProj(), conf.fact[0]->getK(), FULL);
      lat->setPower2(coeffs);
      return lat;
    }

  template<typename Int, typename Real>
    MWCLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf) {
      // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
      // vector format
      NTL::vector<Int> A;
      A.SetLength(conf.fact[0]->getK()+1);
      NTL::clear(A);
      int delay = 0;
      // The program will not run the maxPeriod function if it is not wanted with
      // this condition
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.fact[0]->getK(); i++) A[i+1] = conf.coeff[0][i] * randInt(Int(0), conf.fact[0]->m_MWCb);
        delay++;
        A[0] = -1;
      } while ((A[conf.fact[0]->getK()] == 0));
      return new MWCLattice<Int, Real>(conf.fact[0]->m_MWCb, A, conf.fact[0]->getK(), conf.proj->numProj());
    }

  template<typename Int, typename Real>
    MWCLattice<Int, Real>* nextFullGenerator(ConfigSeek<Int, Real>& conf) {
      Int m(0);
      long exp = conf.fact[0]->getE()-1;
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
      Modulus<Int> mod(conf.fact[0]->getE(), m, true);
      return new MWCLattice<Int, Real>(conf.b, mod.next());
    }

  /*
   * Goin' full random for now
   * */
  template<typename Int, typename Real>
    MMRGLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf) {
      NTL::matrix<Int> A;
      A.SetDims(conf.fact[0]->getK(), conf.fact[0]->getK());
      NTL::clear(A);
      int delay = 0;
      do {
        if (delay >= DELAY) {
          if (timer.timeOver(conf.timeLimit)) return NULL;
          else delay = 0;
        }
        for (long i = 0; i<conf.fact[0]->getK()-1; i++) {
          A[i][i+1] = Int(1);
        }
        for (long i = 0; i<conf.fact[0]->getK(); i++) {
          A[conf.fact[0]->getK()-1][i] = randInt(Int(0), conf.fact[0]->getM());
        }
        delay++;
      } while ((NTL::determinant(A) == 0) || (conf.period[0] && !conf.fact[0]->maxPeriod(A)));
      // Correcting the matrix to a full matrix
      //for (int i = 0; i<conf.fact[0]->getK(); i++) A *= A;
      for (int i = 0; i<conf.fact[0]->getK(); i++)
        for (int j = 0; j<conf.fact[0]->getK(); j++)
          A[i][j] = A[i][j]%conf.fact[0]->getM();
      return new MMRGLattice<Int, Real>(conf.fact[0]->getM(), A, conf.proj->numProj(), conf.fact[0]->getK());
    }

  /**
   * Sing there can be a great variety of combined lattices, this class helps
   * choose the correct constructor for such lattices.
   * */
  template<typename Int, typename Real>
    class ComboLatticeFinder {
      typedef NTL::vector<Int> IntVec;

      public:
      /**
       * The only public and visible function in this class. Calls the correct
       * underlying function.
       * */
      static ComboLattice<Int,Real>* getFunction(ConfigSeek<Int, Real>& conf) {
        return nextGenerator(conf);
      }

      private:

      static ComboLattice<Int, Real>* nextGenerator(ConfigSeek<Int, Real>& conf) {
        // Setting up two vectors. MRGPeriod and MRGLattice do not use the same
        // vector format
        std::vector<std::string> compstr;
        compstr.resize(conf.num_comp);
        for (int k = 0; k < conf.num_comp; k++) {
          // Generating all components
          IntVec A;
          A.SetLength(conf.fact[k]->getK()+1);
          NTL::clear(A);
          int delay = 0;
          std::vector<IntVec> coeffs;
          coeffs.resize(conf.fact[k]->getK());
          // The program will not run the maxPeriod function if it is not wanted with
          // this condition
          do {
            if (delay >= DELAY) {
              if (timer.timeOver(conf.timeLimit)) return NULL;
              else delay = 0;
            }
            for (long i = 0; i < conf.fact[k]->getK(); i++) {
              if (conf.fact[k]->get_type() == MRG) {
                if (conf.search_mode[k] == "pow2") {
                  coeffs[i].resize(NTL::conv<int>(conf.coeff[k][conf.fact[k]->getK()]));
                  A[i+1] = 0;
                  if (conf.coeff[k][i] >= 0) {
                    for (long j = 0; j < conf.coeff[k][conf.fact[k]->getK()]; j++) {
                      coeffs[i][j] = randInt(Int(LatticeTester::Lg(abs(conf.fact[k]->getR()))+1), conf.coeff[k][i]);
                      Int tmp2;
                      NTL::power2(tmp2, coeffs[i][j]);
                      A[i+1] += Int(randInt(0,1)?-1:1) * tmp2;
                    }
                  }
                  auto lat = MRGLattice<Int, Real>(conf.fact[k]->getM(), A, conf.proj->numProj(), conf.fact[k]->getK(), FULL);
                  lat.setPower2(coeffs);
                  compstr[k] = lat.toString();
                } else if (conf.search_mode[k] == "af") {
                  A[i+1] = conf.coeff[k][i] * randInt(Int(0), NTL::SqrRoot(conf.fact[k]->getM()));
                  auto lat = MRGLattice<Int, Real>(conf.fact[k]->getM(), A, conf.proj->numProj(), conf.fact[k]->getK(), FULL);
                  compstr[k] = lat.toString();
                } else if (conf.search_mode[k] == "random") {
                  A[i+1] = randInt(Int(0), conf.coeff[k][i]);
                  auto lat = MRGLattice<Int, Real>(conf.fact[k]->getM(), A, conf.proj->numProj(), conf.fact[k]->getK(), FULL);
                  compstr[k] = lat.toString();
                }
              }
              if (conf.fact[k]->get_type() == MWC)
                A[i+1] = conf.coeff[k][i] * randInt(Int(0), conf.fact[k]->m_MWCb);
            }
            if (conf.fact[k]->get_type() == MWC) {
              A[0] = -1;
              auto lat = MWCLattice<Int, Real>(conf.fact[k]->m_MWCb, A, conf.fact[k]->getK(), conf.proj->numProj());
              compstr[k] = lat.toString();
            }
            delay++;
          } while ((A[conf.fact[k]->getK()] == 0) || (conf.period[k] && !conf.fact[k]->maxPeriod(A)));
          if (conf.fact[k]->get_type() == MRG) {
            IntVec B;
            B.SetLength(conf.fact[k]->getK());
            for (int i = 0; i < conf.fact[k]->getK(); i++) B[i] = A[i+1];
            conf.fact[k]->setA(B);
          } else {
            conf.fact[k]->setA(A);
          }
        }
        MRGLattice<Int, Real>* mrg_lat = getLatCombo<Int, Real>(conf.fact, conf.proj->numProj());
        ComboLattice<Int, Real>* new_lat = new ComboLattice<Int, Real>(conf.fact, *mrg_lat);
        delete mrg_lat;
        for (int i = 0; i < conf.num_comp; i++) new_lat->getCompString(i) = compstr[i];
        return new_lat;
      }
    };
}

template<typename Lat> struct SeekMain {
  typedef typename Lat::Int Int;
  typedef typename Lat::Real Real;
  typedef NTL::vector<Int> IntVec;
  typedef NTL::matrix<Int> IntMat;

  ConfigSeek<Int, Real> conf;
  Lat* lat = 0;

  SeekMain(ConfigSeek<Int, Real>& conf) { this->conf = conf;}

  ~SeekMain() {
    if (lat != NULL) delete lat;
    for (int i = 0; i < conf.num_comp; i++) {
      if (conf.fact[i]) delete conf.fact[i];
    }
    delete conf.proj;
  }

  /*
   * These next function add the tested lattices to the list of the best ones.
   * This only add the lattices that are good enough.
   * */
  void printResults(MeritList<Lat>& bestLattice) {
    *out << "\nSeek: A search program for Random Number Generators\n";
    *out << delim;
    *out << ((conf.num_comp>1)?"Combined generators":"Simple generator")
      << " configuration" << ((conf.num_comp>1)?"s":"") << "\n\n";
    for (int k = 0; k < conf.num_comp; k++) {
      if (k > 0) *out << "\n";
      if (conf.num_comp >1) *out << "Component " << k+1 << ":\n";
      *out << "Generator type: " << toStringGen(conf.fact[k]->get_type()) << "\n";
      if (conf.fact[k]->get_type() == MRG) {
        *out << "Modulo:         m = " << conf.fact[k]->getM() << " = " << conf.fact[k]->getB() << "^"
          << conf.fact[k]->getE();
        if (conf.fact[k]->getR() > 0) *out << "+" << conf.fact[k]->getR();
        if (conf.fact[k]->getR() < 0) *out << conf.fact[k]->getR();
        *out << "\n";
        *out << "Order:          k = " << conf.fact[k]->getK() << "\n";
      } else if (conf.fact[k]->get_type() == MWC) {
      } else if (conf.fact[k]->get_type() == MMRG) {
      }
      *out << (conf.period[0]?"Full":"Any") << " period length\n";
    }
    *out << "\nTest:\n" << (conf.best?"Best":"Worst") << " generators "
      "ranked by ";
    if(conf.criterion == LatticeTester::SPECTRAL) *out
      << (conf.normaType==LatticeTester::NONE?"minimal inverse shortest "
          "non-zero vector length":"Spectral Test") << "\n";
    else if (conf.criterion == LatticeTester::LENGTH) *out << "minimal"
      << " shortest non-zero vector length\n";
    else if (conf.criterion == LatticeTester::BEYER) *out << "their Beyer quotient\n";
    if (conf.normaType != LatticeTester::NONE) {
      *out << "Normalizer used: "
        << LatticeTester::toStringNorma(conf.normaType) << "\n";
    }
    *out << "\nDimensions and projections:\n";
    *out << conf.proj->toString();
    *out << delim;
    if (print_time) {
      *out << "Allowed running time: " << conf.timeLimit << "s.\n";
      *out << "Actual CPU time: " << timer.toString() << "\n";
    }
    *out << "Number of generators kept: " << conf.max_gen << "\n";
    *out << "Number of generators tested: " << conf.num_gen << "\n\n";
    *out << "Retained generators (from best to worst):\n";
    for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
      *out << delim;
      *out << (*it).getLattice() << "\n";
      *out << (*it).toStringMerit() << "\n";
    }
  }

  int print_progress(int old) {
    int per_80 = timer.val(Chrono::SEC)/conf.timeLimit * 80;
    if (per_80 > 80) per_80 = 80;
    if (per_80 < 0) per_80 = 0;
    // We do not print for no reason as this slows the program a lot.
    if (per_80 <= old) return old;
    std::cout << "Program progress: [";
    for (int i = 0; i < per_80; i++) std::cout << "#";
    for (int i = per_80; i < 80; i++) std::cout << " ";
    std::cout << "] ";
    std::cout << std::setw(3) << int(per_80/80.0*100) << " %\r" << std::flush;
    return per_80;
  }

  int Seek (Lat* (*nextGenerator)(ConfigSeek<Int, Real>&) )
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
    // Dynamically allocated objects
    timer.init();

    int old = 0;
    // Launching the tests
    if (conf.progress) {
      old = print_progress(-1);
    }
    MeritList<Lat> bestLattice(conf.max_gen, conf.best);
    do {
      if (lat != NULL) delete lat;
      lat = nextGenerator(conf);
      if (lat == NULL) continue;
      bestLattice.add(test_seek(*lat, conf));
      conf.num_gen++;
      conf.currentMerit = bestLattice.getMerit();
      if (conf.progress) old = print_progress(old);
    } while (!timer.timeOver(conf.timeLimit) && lat);
    std::cout << "\r                                                                                                          \r";
    printResults(bestLattice);
    return 0;
  }
}; // end struct SeekMain
#endif

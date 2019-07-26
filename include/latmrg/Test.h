/*
 * This file contains
 * */
#ifndef LATMRG_TEST_H
#define LATMRG_TEST_H

#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"

#include "latmrg/FigureOfMerit.h"
#include "latmrg/Const.h"

namespace LatMRG {

  /**
   * The methods in this namespace can compute the merit of a reduced lattice
   * for a given projection and return it as a `Dbl`.
   * \todo add extern definitions
   * */
  namespace Merit {

    /**
     * Computes the length of the shortest vector in the lattice.
     * */
    template<typename Lat>
      typename Lat::Float meritL(Lat& lat) {
        typename Lat::IntVec shortest(lat.getBasis()[0]);
        typename Lat::Float tmp;
        LatticeTester::ProdScal<typename Lat::Integ>(shortest, shortest, shortest.length(), tmp);
        return NTL::sqrt(tmp);
      }

    /**
     * Computes the value of the spectral test normalized with `norma`.
     * */
    template<typename Lat>
      typename Lat::Float meritS(Lat& lat,
          LatticeTester::Normalizer<typename Lat::Float>* norma) {
        typename Lat::IntVec shortest(lat.getBasis()[0]);
        typename Lat::Float tmp;
        LatticeTester::ProdScal<typename Lat::Integ>(shortest, shortest, shortest.length(), tmp);
        tmp = NTL::sqrt(tmp)/norma->getBound(shortest.length());
        if (tmp > 1) tmp = typename Lat::Float(1)/tmp;
        return tmp;
      }

    /**
     * Computes the Beyer ratio.
     * \todo implement this
     * */
    template<typename Lat>
      typename Lat::Float meritB(Lat& lat) {
        return typename Lat::Float(0.0);
      }
  }

  /**
   * Functions in this namespace perform lattice reductions.
   * \todo have more flexibility in the reduction parameters so that LLL and BKZ
   * floating point numbers precision can be changed if needed.
   * */
  namespace Reductions {

    /**
     * This performs a BKZ reduction and a shortest vector search.
     * */
    template<typename Int, typename Dbl>
      void reduceFull(LatticeTester::IntLattice<Int, Int, Dbl, Dbl>& lat) {
        LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
        LatticeTester::Reducer<Int, Int, Dbl, Dbl>::maxNodesBB = 100000000000;
        //red.redDieter(0);
        //red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
        red.redLLLNTL(0.999999, LatticeTester::EXPONENT, lat.getDim());
        red.shortestVector(lat.getNorm());
      }

    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);

    /**
     * This performs the BKZ reduction.
     * */
    template<typename Int, typename Dbl>
      void reduceBKZ(LatticeTester::IntLattice<Int, Int, Dbl, Dbl>& lat) {
        LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
        red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
      }

    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);

    /**
     * This performs the LLL reduction.
     * */
    template<typename Int, typename Dbl>
      void reduceLLL(LatticeTester::IntLattice<Int, Int, Dbl, Dbl>& lat) {
        LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
        red.redLLLNTL(0.999999, LatticeTester::EXPONENT, lat.getDim());
      }

    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);

    /**
     * This performs the BKZ reduction before launching a search for a Minkowski
     * reduced basis.
     * */
    template<typename Int, typename Dbl>
      void reduceMink(LatticeTester::IntLattice<Int, Int, Dbl, Dbl>& lat) {
        LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
        red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
        red.reductMinkowski(lat.getDim());
      }

    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLattice<std::int64_t, std::int64_t, double, double>& lat);
    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, double, double>& lat);
    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLattice<NTL::ZZ, NTL::ZZ, NTL::RR, NTL::RR>& lat);
  }

#ifdef LATMRG_USE_CONFIG
  /**
   * This stores the configuration of a problem. This contains many parameters
   * used in the executables and when using the `test()` function.
   *
   * This structure is intended to be included in executables only and needs to
   * be compiled on demand. To do so, simply `#define LATMRG_USE_CONFIG` before
   * including this header.
   * */
  struct Config {
    // For period tests
    DecompType decompm1 = DECOMP, decompr = DECOMP;
    std::string filem1, filer;

    // Data file read parameters
    GenType type;
    // Type of figure of merit
    LatticeTester::NormaType normaType = LatticeTester::NONE;
    LatticeTester::CriterionType criterion;
    LatticeTester::PreReductionType reduction;
    bool use_dual;
#ifdef LATMRG_SEEK
    bool best;
    int num_gen = 0;
    Dbl currentMerit;
#endif
    // Projection
    Projections* proj;

    double timeLimit;
    int max_gen;

    // MRG Specific parameters
    IntVec mult; // MRG multipliers

    // MWC Specific parameters
    Int b; // modulo of MWC recurence

    // MMRG Specific parameters
    IntMat matrix;

    // Combo specific parameters
#ifdef LATMRG_SEEK
    int num_comp;
    std::vector<std::int64_t> comb_order; // k for MRG and MMRG
    IntVec comb_modulo; // m for MRG and MMRG
    // modulo is basis^exponent+rest
    std::vector<std::int64_t> comb_basis;
    std::vector<std::int64_t> comb_exponent;
    std::vector<std::int64_t> comb_rest;
    std::vector<bool> comb_period;
    std::vector<MRGComponent<Int>*> comb_fact;
#endif

    // Shared components names
    std::int64_t order; // k for MRG and MMRG
    Int modulo; // m for MRG and MMRG
    // modulo is basis^exponent+rest
    std::int64_t basis;
    std::int64_t exponent;
    std::int64_t rest;
    bool period;
  };

#ifdef LATMRG_USE_TEST
  /**
   * This performs a test on a lattice. That is, given a lattice and a `Config`
   * it populates a `FigureOfMerit` objects and returns it.
   * */
  template<typename Lat>
    FigureOfMerit<Lat> test(Lat & lattice, Config& conf) {
      typedef typename Lat::Integ Int;
      typedef typename Lat::IntVec IntVec;
      typedef typename Lat::Float Dbl;
      Projections* proj(conf.proj);

      LatticeTester::Normalizer<Dbl>* norma = lattice.getNormalizer(conf.normaType, 0, true);
      std::vector<Dbl> results;
      std::vector<IntVec> vectors;
      lattice.buildBasis(proj->minDim());
      for (int i = proj->minDim(); i <= proj->maxDim(); i++){
#ifdef LATMRG_LAT
        if (!((i-1)%5)) std::cout << "i " << i-1 << " time " << timer.val(LatMRG::Chrono::SEC) << "\n";
        if (timer.timeOver(conf.timeLimit)) {
          std::cout << "On projection " << i << std::endl;
          FigureOfMerit<Lat> figure(lattice, *proj);
          figure.addMerit(results, vectors);
          figure.setFinished();
          figure.computeMerit("min");

          delete norma;
          return figure;
        }
#endif
        // Changing to the dual
        if (conf.use_dual) lattice.dualize();
        // Reducing the lattice
        if (conf.reduction == LatticeTester::FULL)
          Reductions::reduceFull(lattice);
        else if (conf.reduction == LatticeTester::LLL)
          Reductions::reduceLLL(lattice);
        else if (conf.reduction == LatticeTester::BKZ)
          Reductions::reduceBKZ(lattice);
        else if (conf.reduction == LatticeTester::NOPRERED)
          Reductions::reduceMink(lattice);
        // Computing the merit of the lattice
        Dbl tmp;
        if (conf.criterion == LatticeTester::LENGTH) tmp = Merit::meritL(lattice);
        if (conf.criterion == LatticeTester::SPECTRAL) tmp = Merit::meritS(lattice, norma);
        if (conf.criterion == LatticeTester::BEYER) tmp = Merit::meritB(lattice);
        results.push_back(tmp);
        vectors.push_back(lattice.getBasis()[0]);
#ifdef LATMRG_SEEK
        // Rejecting lattices that won't make it
        if (tmp < conf.currentMerit) {
          delete norma;
          return FigureOfMerit<Lat>(lattice, *conf.proj);
        }
#endif
        // Changing back to the primal and increasing the dimension
        if (conf.use_dual) lattice.dualize();
        if (proj->minDim() < proj->maxDim()) lattice.incDim();
      }
#ifdef LATMRG_LAT
        std::cout << "Seq time " << timer.val(LatMRG::Chrono::SEC) << "\n";
#endif

      // Testing projections if there are anyo
      // This is done separately because sequential testing is much more efficient
      for (int i = 2; i <= proj->numProj(); i++) {
        proj->resetDim(i);
        lattice.buildBasis(proj->projDim()[i-1]+1);
        while(!proj->end(1)) {
          // Building the projection
          LatticeTester::IntLattice<Int, Int, Dbl, Dbl> proj_lat(lattice.getModulo(), lattice.getOrder(), i, true);
          LatticeTester::Coordinates iter(proj->next());
#ifdef LATMRG_LAT
          if (timer.timeOver(conf.timeLimit)) {
            std::cout << "On projection " << iter << std::endl;
            FigureOfMerit<Lat> figure(lattice, *proj);
            figure.addMerit(results, vectors);
            figure.setFinished();
            figure.computeMerit("min");

            delete norma;
            return figure;
          }
#endif
          lattice.buildProjection(&proj_lat, iter);
          norma->setLogDensity(Dbl(-i*log(lattice.getModulo())
                +log(abs(NTL::determinant(proj_lat.getBasis())))));
          if (conf.use_dual) proj_lat.dualize();
          // Reduction
          if (conf.reduction == LatticeTester::FULL)
            Reductions::reduceFull(proj_lat);
          else if (conf.reduction == LatticeTester::LLL)
            Reductions::reduceLLL(proj_lat);
          else if (conf.reduction == LatticeTester::BKZ)
            Reductions::reduceBKZ(proj_lat);
          else if (conf.reduction == LatticeTester::NOPRERED)
            Reductions::reduceMink(proj_lat);

          // Figure of merit
          Dbl tmp;
          if (conf.criterion == LatticeTester::LENGTH) tmp = Merit::meritL(proj_lat);
          else if (conf.criterion == LatticeTester::SPECTRAL) tmp = Merit::meritS(proj_lat, norma);
          else if (conf.criterion == LatticeTester::BEYER) tmp = Merit::meritB(proj_lat);
          results.push_back(tmp);
          vectors.push_back(proj_lat.getBasis()[0]);
#ifdef LATMRG_SEEK
          // Rejecting lattices that won't make it
          if (tmp < conf.currentMerit) {
            delete norma;
            return FigureOfMerit<Lat>(lattice, *proj);
          }
#endif
        }
#ifdef LATMRG_LAT
        std::cout << "dim " << i << " time " << timer.val(LatMRG::Chrono::SEC) << "\n";
#endif
      }

      FigureOfMerit<Lat> figure(lattice, *proj);
      figure.addMerit(results, vectors);
      figure.setFinished();
      figure.computeMerit("min");

      delete norma;
      return figure;
    }
#endif
#endif

} // end namespace LatMRG
#endif

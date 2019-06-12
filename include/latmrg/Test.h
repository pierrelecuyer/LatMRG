#ifndef LATMRG_TEST_H
#define LATMRG_TEST_H

namespace LatMRG {

  /**
   * The next few methods compute a figure of merit depending on the specified
   * figure of merit. The lattice has to be properly reduced when calling this.
   * */
  // Shortest vector length
  Dbl meritL(IntLattice<Int, Int, Dbl, Dbl>& lat,
      Normalizer<Dbl>* norma) {
    IntVec shortest(lat.getBasis()[0]);
    Dbl tmp;
    LatticeTester::ProdScal<Int>(shortest, shortest, shortest.length(), tmp);
    return NTL::sqrt(tmp);
  }

  // Spectral test
  Dbl meritS(IntLattice<Int, Int, Dbl, Dbl>& lat,
      Normalizer<Dbl>* norma) {
    IntVec shortest(lat.getBasis()[0]);
    Dbl tmp;
    LatticeTester::ProdScal<Int>(shortest, shortest, shortest.length(), tmp);
    tmp = NTL::sqrt(tmp)/norma->getBound(shortest.length());
    if (tmp > 1) tmp = Dbl(1)/tmp;
    return tmp;
  }

  // Beyer ratio
  Dbl meritB(IntLattice<Int, Int, Dbl, Dbl>& lat,
      Normalizer<Dbl>* norma) {
    return 0;
  }

  void reduceFull(IntLattice<Int, Int, Dbl, Dbl>& lat) {
    LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
    red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
    red.shortestVector(lat.getNorm());
  }

  void reduceBKZ(IntLattice<Int, Int, Dbl, Dbl>& lat) {
    LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
    red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lat.getDim());
  }

  void reduceLLL(IntLattice<Int, Int, Dbl, Dbl>& lat) {
    LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
    red.redLLLNTL(0.999999, LatticeTester::QUADRUPLE, lat.getDim());
  }

  void reduceMink(IntLattice<Int, Int, Dbl, Dbl>& lat) {
    LatticeTester::Reducer<Int, Int, Dbl, Dbl> red(lat);
    red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, lat.getDim());
    red.reductMinkowski(lat.getDim());
  }

  // The configuration of a problem
  struct Config {
    // For period tests
    DecompType decompm1 = DECOMP, decompr = DECOMP;
    std::string filem1, filer;
    TestList* bestLattice;

    // Data file read parameters
    GenType type;
    // Type of figure of merit
    LatticeTester::NormaType normaType = LatticeTester::NONE;
    LatticeTester::CriterionType criterion;
    LatticeTester::PreReductionType reduction;
    bool use_dual;
#ifdef LATMRG_SEEK
    bool best;
#endif
    // Projection
    int numProj;
    int minDim, maxDim;
    std::vector<std::size_t> projDim;
    Projections* proj;

    double timeLimit;
    int num_gen;

    // MRG Specific parameters
    IntVec mult; // MRG multipliers

    // MWC Specific parameters
    Int b; // modulo of MWC recurence

    // Shared components names
    std::int64_t order; // k for MRG and MMRG
    Int modulo; // m for MRG and MMRG
    // modulo is basis^exponent+rest
    std::int64_t basis;
    std::int64_t exponent;
    std::int64_t rest;
    bool period;
  };

  DblVec test(IntLattice<Int, Int, Dbl, Dbl> & lattice, Config& conf) {
    LatticeTester::Normalizer<Dbl>* norma = lattice.getNormalizer(conf.normaType, 0, true);
    DblVec results(0);
    lattice.buildBasis(conf.minDim);
    for (int i = conf.minDim; i <= conf.maxDim; i++){
      // Changing to the dual
      if (conf.use_dual) lattice.dualize();
      // Reducing the lattice
      if (conf.reduction == LatticeTester::FULL)
        reduceFull(lattice);
      else if (conf.reduction == LatticeTester::LLL)
        reduceLLL(lattice);
      else if (conf.reduction == LatticeTester::BKZ)
        reduceBKZ(lattice);
      else if (conf.reduction == LatticeTester::NOPRERED)
        reduceMink(lattice);
      // Computing the merit of the lattice
      Dbl tmp;
      if (conf.criterion == LatticeTester::LENGTH) tmp = meritL(lattice, norma);
      if (conf.criterion == LatticeTester::SPECTRAL) tmp = meritS(lattice, norma);
      if (conf.criterion == LatticeTester::BEYER) tmp = meritB(lattice, norma);
      results.append(tmp);
#ifdef LATMRG_SEEK
      // Rejecting lattices that won't make it
      if (tmp < conf.bestLattice->getMerit()) {
        results[0] = 1-conf.best;
        return results;
      }
#endif
      // Changing back to the primal and increasing the dimension
      if (conf.use_dual) lattice.dualize();
      lattice.incDim();
    }

    // Testing projections if there are any
    for (int i = 2; i <= conf.numProj; i++) {
      conf.proj->resetDim(i);
      lattice.buildBasis(conf.projDim[i-1]+1);
      while(!conf.proj->end(1)) {
        // Building the projection
        IntLattice<Int, Int, Dbl, Dbl> proj_lat(conf.modulo, conf.order, i, true);
        LatticeTester::Coordinates iter(conf.proj->next());
        lattice.buildProjection(&proj_lat, iter);
        norma->setLogDensity(Dbl(-i*log(conf.modulo)
              +log(abs(NTL::determinant(proj_lat.getBasis())))));
        if (conf.use_dual) proj_lat.dualize();
        // Reduction
        if (conf.reduction == LatticeTester::FULL)
          reduceFull(proj_lat);
        else if (conf.reduction == LatticeTester::LLL)
          reduceLLL(proj_lat);
        else if (conf.reduction == LatticeTester::BKZ)
          reduceBKZ(proj_lat);
        else if (conf.reduction == LatticeTester::NOPRERED)
          reduceMink(proj_lat);
        // Figure of merit
        Dbl tmp;
        if (conf.criterion == LatticeTester::LENGTH) tmp = meritL(proj_lat, norma);
        else if (conf.criterion == LatticeTester::SPECTRAL) tmp = meritS(proj_lat, norma);
        else if (conf.criterion == LatticeTester::BEYER) tmp = meritB(proj_lat, norma);
        results.append(tmp);
#ifdef LATMRG_SEEK
        // Rejecting lattices that won't make it
        if (tmp < conf.bestLattice->getMerit()) {
          results[0] = 1-conf.best;
          return results;
        }
#endif
      }
    }

    delete norma;
    return results;
  }

} // end namespace LatMRG
#endif

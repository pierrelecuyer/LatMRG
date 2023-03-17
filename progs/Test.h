typedef typename Lat::Int Int;
typedef typename Lat::IntVec IntVec;
typedef typename Lat::Real Real;
Projections* proj(conf.proj);

LatticeTester::Normalizer<Real>* norma = lattice.getNormalizer(conf.normaType, 0, conf.use_dual);
std::vector<Real> results;
std::vector<IntVec> vectors;
lattice.buildBasis(proj->minDim());
for (int i = proj->minDim(); i <= proj->maxDim(); i++){
#ifdef LATMRG_LAT
  // if (!((i-1)%5)) std::cout << "i " << i-1 << " time " << timer.val(LatMRG::Chrono::SEC) << "\n";
  // if (timer.timeOver(conf.timeLimit)) {
  //   std::cout << "On projection " << i << std::endl;
  //   FigureOfMerit<Lat> figure(lattice, *proj);
  //   figure.addMerit(results, vectors);
  //   figure.setFinished();
  //   figure.computeMerit("min");

  //   delete norma;
  //   return figure;
  // }
#endif
  // Changing to the dual
  if (conf.use_dual) lattice.dualize();
  // Reducing the lattice
  if (conf.reduction == LatticeTester::FULL)
    Reductions::reduceFull(lattice, conf.nodesBB);
  else if (conf.reduction == LatticeTester::LLL)
    Reductions::reduceLLL(lattice);
  else if (conf.reduction == LatticeTester::BKZ)
    Reductions::reduceBKZ(lattice);
  else if (conf.reduction == LatticeTester::NOPRERED)
    Reductions::reduceMink(lattice);
  // Computing the merit of the lattice
  Real tmp;
  if (conf.criterion == LatticeTester::LENGTH) tmp = Merit::meritL(lattice);
  if (conf.criterion == LatticeTester::SPECTRAL) tmp = Merit::meritS(lattice, norma);
  if (conf.criterion == LatticeTester::BEYER) tmp = Merit::meritB(lattice);
  results.push_back(tmp);
  vectors.push_back(lattice.getBasis()[0]);
#ifdef LATMRG_SEEK
  // Rejecting lattices that won't make it
  if ((tmp < conf.currentMerit) && conf.best) {
    delete norma;
    return FigureOfMerit<Lat>(lattice, *conf.proj);
  }
#endif
  // Changing back to the primal and increasing the dimension
  if (conf.use_dual) lattice.dualize();
  if (proj->minDim() < proj->maxDim()) {
    lattice.incDim();
  }
}
#ifdef LATMRG_LAT
// std::cout << "Seq time " << timer.val(LatMRG::Chrono::SEC) << "\n";
#endif

// Testing projections if there are anyo
// This is done separately because sequential testing is much more efficient
for (int i = 2; i <= proj->numProj(); i++) {
  proj->resetDim(i);
  lattice.buildBasis(proj->projDim()[i-1]+1);
  while(!proj->end(1)) {
    // Building the projection
    LatticeTester::IntLatticeExt<Int, Real> proj_lat(lattice.getModulo(), lattice.getOrder(), i, true);
    LatticeTester::Coordinates iter(proj->next());
#ifdef LATMRG_LAT
    // if (timer.timeOver(conf.timeLimit)) {
    //   std::cout << "On projection " << iter << std::endl;
    //   FigureOfMerit<Lat> figure(lattice, *proj);
    //   figure.addMerit(results, vectors);
    //   figure.setFinished();
    //   figure.computeMerit("min");

    //   delete norma;
    //   return figure;
    // }
#endif
    lattice.buildProjection(&proj_lat, iter);
    norma->setLogDensity(Real(-i*log(lattice.getModulo())
          +log(abs(NTL::determinant(proj_lat.getBasis())))));
    if (conf.use_dual) proj_lat.dualize();
    // Reduction
    if (conf.reduction == LatticeTester::FULL)
      Reductions::reduceFull(proj_lat, conf.nodesBB);
    else if (conf.reduction == LatticeTester::LLL)
      Reductions::reduceLLL(proj_lat);
    else if (conf.reduction == LatticeTester::BKZ)
      Reductions::reduceBKZ(proj_lat);
    else if (conf.reduction == LatticeTester::NOPRERED)
      Reductions::reduceMink(proj_lat);

    // Figure of merit
    Real tmp;
    if (conf.criterion == LatticeTester::LENGTH) tmp = Merit::meritL(proj_lat);
    else if (conf.criterion == LatticeTester::SPECTRAL) tmp = Merit::meritS(proj_lat, norma);
    else if (conf.criterion == LatticeTester::BEYER) tmp = Merit::meritB(proj_lat);
    results.push_back(tmp);
    vectors.push_back(proj_lat.getBasis()[0]);
#ifdef LATMRG_SEEK
    // Rejecting lattices that won't make it
    if ((tmp < conf.currentMerit) && conf.best) {
      delete norma;
      return FigureOfMerit<Lat>(lattice, *proj);
    }
#endif
  }
#ifdef LATMRG_LAT
  // std::cout << "dim " << i << " time " << timer.val(LatMRG::Chrono::SEC) << "\n";
#endif
}

FigureOfMerit<Lat> figure(lattice, *proj);
figure.addMerit(results, vectors);
figure.setFinished();
figure.computeMerit("min");

delete norma;
return figure;

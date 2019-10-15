/*
 * This file contains implementations and configurations for tests.
 * The main functions of this file are only defined when certain preprocessing
 * macros exist.
 * */
#ifndef LATMRG_TEST_H
#define LATMRG_TEST_H

#include "latticetester/IntLattice.h"
#include "latticetester/Reducer.h"

#include "latmrg/FigureOfMerit.h"
#include "latmrg/Const.h"
#include "latmrg/MRGComponent.h"

namespace LatMRG {

  /**
   * The methods in this namespace can compute the merit of a reduced lattice
   * for a given projection and return it as a `Dbl`.
   * */
  namespace Merit {

    /**
     * Computes the length of the shortest vector in the lattice.
     * */
    template<typename Lat>
      typename Lat::Dbl meritL(Lat& lat) {
        typename Lat::IntVec shortest(lat.getBasis()[0]);
        typename Lat::Dbl tmp;
        LatticeTester::ProdScal<typename Lat::Int>(shortest, shortest, shortest.length(), tmp);
        return NTL::sqrt(tmp);
      }

    /**
     * Computes the value of the spectral test normalized with `norma`.
     * */
    template<typename Lat>
      typename Lat::Dbl meritS(Lat& lat,
          LatticeTester::Normalizer<typename Lat::Dbl>* norma) {
        typename Lat::IntVec shortest(lat.getBasis()[0]);
        typename Lat::Dbl tmp;
        LatticeTester::ProdScal<typename Lat::Int>(shortest, shortest, shortest.length(), tmp);
        tmp = NTL::sqrt(tmp)/norma->getBound(shortest.length());
        if (tmp > 1) tmp = typename Lat::Dbl(1)/tmp;
        return tmp;
      }

    /**
     * Computes the Beyer ratio.
     * \todo implement this
     * */
    template<typename Lat>
      typename Lat::Dbl meritB(Lat& lat) {
        return typename Lat::Dbl(0.0);
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

} // end namespace LatMRG
#endif

/*
 * This file contains implementations and configurations for tests.
 * The main functions of this file are only defined when certain preprocessing
 * macros exist.
 * */
#ifndef LATMRG_TEST_H
#define LATMRG_TEST_H

#include "latticetester/IntLatticeExt.h"
#include "latticetester/Reducer.h"

#include "latmrg/FigureOfMerit.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/MRGPeriod.h"

namespace LatMRG {

  /**
   * The methods in this namespace can compute the merit of a reduced lattice
   * for a given projection and return it as a `Real`.
   * */
  namespace Merit {

    /**
     * Computes the length of the shortest vector in the lattice.
     * */
    template<typename Lat>
      typename Lat::Real meritL(Lat& lat) {
        // typename Lat::IntVec shortest(lat.getBasis()[0]);
        // typename Lat::Real tmp;
        // LatticeTester::ProdScal<typename Lat::Int>(shortest, shortest, shortest.length(), tmp);
        lat.updateVecNorm();
        if (lat.getNorm() == LatticeTester::L2NORM) return NTL::sqrt(lat.getVecNorm(0));
        else return lat.getVecNorm(0);
      }

    /**
     * Computes the value of the spectral test normalized with `norma`.
     * */
    template<typename Lat>
      typename Lat::Real meritS(Lat& lat,
          LatticeTester::Normalizer<typename Lat::Real>* norma) {
        typename Lat::IntVec shortest(lat.getBasis()[0]);
        typename Lat::Real tmp;
        LatticeTester::ProdScal<typename Lat::Int>(shortest, shortest, shortest.length(), tmp);
        tmp = NTL::sqrt(tmp)/norma->getBound(shortest.length());
        if (tmp > 1) tmp = typename Lat::Real(1)/tmp;
        return tmp;
      }

    /**
     * Computes the Beyer ratio.
     * \todo implement this
     * */
    template<typename Lat>
      typename Lat::Real meritB(Lat& lat) {
        return typename Lat::Real(0.0);
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
    template<typename Int, typename Real>
      void reduceFull(LatticeTester::IntLatticeExt<Int, Int, Real, Real>& lat, std::int64_t maxNodesBB = 100000000000) {
        LatticeTester::Reducer<Int, Int, Real, Real> red(lat);
        LatticeTester::Reducer<Int, Int, Real, Real>::maxNodesBB = maxNodesBB;
        //red.redDieter(0);
        //red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
        red.redLLLNTL(0.999999, LatticeTester::EXPONENT, lat.getDim());
        red.shortestVector(lat.getNorm());
      }

    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLatticeExt<std::int64_t, std::int64_t, double>& lat, std::int64_t);
    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, double>& lat, std::int64_t);
    /**
     * Instanciation of `reduceFull()`.
     * */
    extern template void reduceFull(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, NTL::RR>& lat, std::int64_t);

    /**
     * This performs the BKZ reduction.
     * */
    template<typename Int, typename Real>
      void reduceBKZ(LatticeTester::IntLatticeExt<Int, Int, Real, Real>& lat) {
        LatticeTester::Reducer<Int, Int, Real, Real> red(lat);
        red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
      }

    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLatticeExt<std::int64_t, std::int64_t, double>& lat);
    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, double>& lat);
    /**
     * Instanciation of `reduceBKZ()`.
     * */
    extern template void reduceBKZ(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, NTL::RR>& lat);

    /**
     * This performs the LLL reduction.
     * */
    template<typename Int, typename Real>
      void reduceLLL(LatticeTester::IntLatticeExt<Int, Int, Real, Real>& lat) {
        LatticeTester::Reducer<Int, Int, Real, Real> red(lat);
        red.redLLLNTL(0.999999, LatticeTester::EXPONENT, lat.getDim());
      }

    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLatticeExt<std::int64_t, std::int64_t, double>& lat);
    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, double>& lat);
    /**
     * Instanciation of `reduceLLL()`.
     * */
    extern template void reduceLLL(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, NTL::RR>& lat);

    /**
     * This performs the BKZ reduction before launching a search for a Minkowski
     * reduced basis.
     * */
    template<typename Int, typename Real>
      void reduceMink(LatticeTester::IntLatticeExt<Int, Int, Real, Real>& lat) {
        LatticeTester::Reducer<Int, Int, Real, Real> red(lat);
        red.redBKZ(0.999999, 10, LatticeTester::EXPONENT, lat.getDim());
        red.reductMinkowski(lat.getDim());
      }

    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLatticeExt<std::int64_t, std::int64_t, double>& lat);
    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, double>& lat);
    /**
     * Instanciation of `reduceMink()`.
     * */
    extern template void reduceMink(LatticeTester::IntLatticeExt<NTL::ZZ, NTL::ZZ, NTL::RR>& lat);
  }

} // end namespace LatMRG
#endif

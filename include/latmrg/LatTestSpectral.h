#ifndef LATTESTSPECTRAL_H
#define LATTESTSPECTRAL_H

#include <cmath>

#include "latticetester/Reducer.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Const.h"

#include "latmrg/LatticeTest.h"
#include "latmrg/Merit.h"
#include "latmrg/PolyPE.h"
#include "latmrg/LatticeTest.h"

namespace LatMRG {

  /**
   * This class implements the *spectral* test. It implements the abstract
   * class `LatticeTest`. The figure of merit for this test is the length of
   * the shortest vector in the *primal* or in the *dual* lattice computed with
   * different norms. For the standard spectral test, the figure of merit is
   * based on the length of the shortest non-zero vector in the *dual* lattice,
   * using the \f${\mathcal{L}}_2\f$ norm to compute the length of vectors, and
   * the inverse of this length gives the maximal distance between successive
   * hyperplanes covering all the points in the *primal* lattice. If one
   * computes the length of the shortest non-zero vector in the *dual* lattice
   * using the \f${\mathcal{L}}_1\f$ norm instead, one obtains the minimal
   * number of hyperplanes covering all the points of the *primal* lattice.
   *
   * The main program is obtained by compiling the `LatMain.cc` file (see the
   * description of module `LatMain`, on page (FIXME: page#), for how to run a
   * program).
   *
   */
  /**
   * Currently, when the spectral test is applied in multiple dimensions, the
   * dual is switched with the primal at each step, expanded re switched to get
   * the dual into the reducer and the figure of merit is calculated. There is
   * no check to see if the dual is present even though it is needed. LatTestSpectral
   * asks for an IntLattice even though the incDim method in the IntLattice base
   * class is very primitive. This is kind of a miracle that this works.
   * */
  template<typename Int, typename Dbl>
    class LatTestSpectral : public LatticeTest<Int, Dbl> {
      public:

        /**
         * Constructor. The *spectral* test will be applied to the lattice `lat`
         * using normalizer `normal` to normalize the figure of merit.
         * 
         * `lat` should point to a subclass of IntLattice!!!
         */
        LatTestSpectral (LatticeTester::Normalizer<Dbl> * normal,
            LatticeTester::IntLattice<Int, Int, Dbl, Dbl> * lat);

        /**
         * Destructor.
         */
        ~LatTestSpectral ();

        /**
         * Applies the spectral test for dimensions varying from `fromDim` to
         * `toDim`. Whenever the normalized value of the merit is smaller than
         * `minVal` for any dimension, the method returns `false` immediately.
         * The method returns `false` if the test was interrupted for any
         * reason before completion, and it returns `true` upon success. The
         * results of the last test are kept in <tt>m_merit</tt>.
         */
        bool test (int fromDim, int toDim, double minVal[]);

        /**
         * Similar to `test` above, but with the weights `weights`. `weights`
         * is the array of the weights of all projections defined as follows:
         *
         * `weights[0]` is the weight of projections [1, ... ,
         * <tt>fromDim</tt>];<br><tt>weights[1]</tt> is the weight of
         * projections [1, ... , `fromDim` + 1];<br>
         * \f$\vdots\f$ <br><tt>weights</tt>[<tt>toDim - fromDim</tt>] is the
         * weight of projections [1, ... , <tt>toDim</tt>].
         *
         * If `weights` = 0, it means unit weight for all projections.
         */
        bool test (int fromDim, int toDim, double minVal[], const double* weights);

        /**
         * Performs the spetral test exactly as LatTestSpectral::test, except it
         * does not perform branch and bound. There are currently two options for
         * `speed`. If `speed==1` then the figures of merit are approximated with
         * the BKZ pre-reduction and if `speed==2` the figures of merit are 
         * approximated with the LLL pre-reduction. The point of this method is 
         * to speed up research of good RNGs. RNGs that have okay-ish figures of
         * merit after using only an approximated spectral test can then be tested
         * more extensively on many projections and dimensions.
         * */
        bool quicktest(int fromDim, int toDim, double minVal[], int speed);

        /**
         * Sets the lower bound on the square length of the shortest vector in
         * each dimension, based on the spectral value `S2`.
         */
        void setLowerBoundL2 (double S2);

        /**
         * Similar to `setLowerBoundL2` above, but with the weights `weights`.
         */
        void setLowerBoundL2 (double S2, const double* weights);

        /**
         * Returns the normalizer used in this test.
         */
        const LatticeTester::Normalizer<Dbl>* getNormalizer() { return this->m_normalizer; }

      private:

        /**
         * The normalizer used to normalize the figure of merit.
         */
        LatticeTester::Normalizer<Dbl>* m_normalizer;

        /**
         * The lower bound on the square length of the shortest vector in each
         * dimension. As soon as a vector of length smaller than this bound is
         * found, the search for the shortest vector in this lattice is stopped
         * and the lattice is rejected.
         */
        NTL::vector<Dbl> m_boundL2;

        /**
         * Initializes the constants <tt>m_S2toL2</tt> below, necessary to
         * compute the lower bounds in `setLowerBoundL2`, for all dimensions
         * \f$d\f$ such that `dim1` \f$\le d \le\f$ `dim2`. This function must
         * be called only after the lattice has been built (after a call to
         * `buildBasis` or more specifically `initStates`, since the constants
         * depend on the initialization in <tt>initStates</tt>).
         */
        void initLowerBoundL2 (int dim1, int dim2);

        /**
         * These precomputed constants allows the calculation of the square
         * length of a lattice vector \f$\ell_2\f$ from a value of the merit
         * \f$S_2\f$ for each dimension \f$i\f$, i.e. \f$\ell_2[i] =
         * S_2[i]*\f$<tt>m_S2toL2[</tt>\f$i\f$<tt>]</tt>.
         */
        double *m_S2toL2;

        /**
         * Prepares and dispatches the results for dimension `dim` to all
         * observers attached to this test.
         */
        void prepAndDisp (int dim);

        /**
         * Sends the initialization message to all observers attached to this
         * test.
         */
        void init ();
    };

  //===========================================================================

  template<typename Int, typename Dbl>
    LatTestSpectral<Int, Dbl>::LatTestSpectral (LatticeTester::Normalizer<Dbl> * normal, LatticeTester::IntLattice<Int, Int, Dbl, Dbl> * lat): LatticeTest<Int, Dbl> (lat)
  {
    this->m_criter = LatticeTester::SPECTRAL;
    this->m_normalizer = normal;
  }

  //===========================================================================

  template<typename Int, typename Dbl>
    LatTestSpectral<Int, Dbl>::~LatTestSpectral ()
    {
      this->m_boundL2.kill();
      delete [] this->m_S2toL2;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void LatTestSpectral<Int, Dbl>::initLowerBoundL2 (int dim1, int dim2)
    {
      if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
        for (int i = dim1; i <= dim2; i++) {
          this->m_S2toL2[i] = this->m_normalizer->getGamma (i);
          this->m_S2toL2[i] *= exp2 (this->m_lat->getLgVolDual2 (i) / i);
        }

      } else if (this->m_lat->getNorm () == LatticeTester::L1NORM) {
        if (this->m_dualF) {
          for (int i = dim1; i <= dim2; i++)
            this->m_S2toL2[i] = trunc(exp2 ((this->m_lat->getLgVolDual2 (i) / 2.0
                    + this->m_normalizer->getGamma(i)) / i));
        } else {
          // Je ne suis pas sûr que ce soit correct pour le primal
          for (int i = dim1; i <= dim2; i++)
            this->m_S2toL2[i] = exp2 ((this->m_lat->getLgVolDual2 (i) / 2.0
                  + this->m_normalizer->getGamma(i)) / i);
        }

      } else {
        for (int i = dim1; i <= dim2; i++) {
          this->m_S2toL2[i] = 1.0;
        }
      }
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void LatTestSpectral<Int, Dbl>::setLowerBoundL2 (double S2)
    {
      setLowerBoundL2(S2, 0);
    }

  template<typename Int, typename Dbl>
    void LatTestSpectral<Int, Dbl>::setLowerBoundL2 (double S2, const double* weights)
    {
      int i;
      Dbl m2;
      NTL::conv(m2, this->m_lat->getModulo());
      m2 = m2 * m2;
      if (S2 <= 0.0) {
        for (i = this->m_fromDim; i <= this->m_toDim; i++)
          this->m_boundL2[i] = 0;
        return;
      }

      if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
        if (!this->m_dualF) {
          // On a multiplié les coordonnées du primal par m pour calculer
          // uniquement avec des entiers: meilleure précision.
          for (i = this->m_fromDim; i <= this->m_toDim; i++)
            NTL::conv(this->m_boundL2[i],m2 * S2 * this->m_S2toL2[i]);

        } else {
          for (i = this->m_fromDim; i <= this->m_toDim; i++)
            NTL::conv(this->m_boundL2[i], S2*this->m_S2toL2[i]);
        }

      } else {
        if (!this->m_dualF) {
          for (i = this->m_fromDim; i <= this->m_toDim; i++) {
            NTL::conv(this->m_boundL2[i], S2*this->m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            this->m_boundL2[i] = this->m_boundL2[i]*this->m_boundL2[i];
            // On a multiplié les coordonnées du primal par m pour calculer
            // uniquement avec des entiers: meilleure précision.
            this->m_boundL2[i] = this->m_boundL2[i] * m2;
          }

        } else {
          for (i = this->m_fromDim; i <= this->m_toDim; i++) {
            NTL::conv(this->m_boundL2[i], S2*this->m_S2toL2[i]);
            // In Reducer, we compare the square length of vectors
            this->m_boundL2[i] = this->m_boundL2[i]*this->m_boundL2[i];
          }
        }
      }

      // Merit values will be scaled by the inverse weight after normalization.
      // It is the scaled values that need to be compared, but the Reducer doesn't know
      // about scaling. So we scale the bounds below by the weight, which ensures the
      // ordering of the scaled values is preserved.
      // FIXME: should it be the squared weight for LatticeTester::L1NORM?
      if (weights) {
        for (i = this->m_fromDim; i <= this->m_toDim; i++)
          this->m_boundL2[i] *= weights[i];
      }
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    bool LatTestSpectral<Int, Dbl>::test (int fromDim, int toDim, double minVal[])
    {
      return test (fromDim, toDim, minVal, 0);
    }

  // weights is the array of the weights of all projections as follows:
  //    weights[i] is the weight for projection [1, ..., i]
  // weights == 0 means unit weight for all projections

  //===========================================================================

  template<typename Int, typename Dbl>
    bool LatTestSpectral<Int, Dbl>::quicktest (int fromDim, int toDim, double minVal[],
        int speed)
    {
      double* weights = NULL;
      //BOOST_DISPLAY

      this->m_merit.setDim(toDim);
      this->m_fromDim = fromDim;
      this->m_toDim = toDim;
      init ();
      this->resetFromDim (this->m_lat->getOrder (), fromDim);
      if (this->m_normalizer->getNorm () != this->m_lat->getNorm ()) {
        std::cout << "SpectralTest: conflict between NormType and Normalizer" <<
          std::endl;
        exit (EXIT_FAILURE);
      }

      double mr;
      NTL::conv (mr, this->m_lat->getModulo ());
      double temp;

      //Dbl te;
      //double lgvv1 = 0.0;
      // "unused variables" according to the compiler

      while (this->m_lat->getDim () < fromDim)
        this->m_lat->incDim ();

      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red (*this->m_lat);

      if (this->m_S2toL2[fromDim] <= 0.0)
        initLowerBoundL2 (fromDim, toDim);
      setLowerBoundL2 (minVal[toDim], weights);   // same S2 for all dim
      red.setBoundL2 (this->m_boundL2, fromDim, toDim);

      while (true) {

        if (this->m_dualF)
          this->m_lat->dualize ();

        int dim = this->m_lat->getDim ();

        // pre-reduction step before BB with default parameters
        if (speed == 1) {
          red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, dim);
        } else if (speed == 2) {
          red.redLLLNTL(0.999999, LatticeTester::QUADRUPLE, dim);
        } else {
          std::cout << "Wrong speed selection\n";
          return false;
        }
        //PW_TODO: hard coding?

        //std::cout << "dual basis = \n" << this->m_lat->getBasis() << std::endl;
        //std::cout << "primal basis = \n" << this->m_lat->getDualBasis() << std::endl;


        // Calcul de D2. Pour Norm # L2NORM, suppose que VV est a jour.
        if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
          this->m_lat->updateScalL2Norm (0);
          NTL::conv (temp, this->m_lat->getVecNorm (0));
        } else {
          NTL::conv (temp, red.getMinLength ());
        }
        if (3 == this->m_detailF) {
          this->dispatchLatUpdate(*this->m_lat);
        }
        // weight factor:
        // require a higher merit for more important projections by dividing
        // their merit by the weight of the projection
        double weight = weights ? weights[dim] : 1.0;

        this->m_merit.getMerit (dim) = temp / weight;

        if (dim <= toDim) { // Calcul de S2.

          //std::cout << "dim = " << dim << std::endl;
          //std::cout << "order = " << this->m_lat->getOrder() << std::endl;
          //std::cout << "density AVANT = " << exp(this->m_normalizer->getLogDensity()) << std::endl;

          //updating value of matrix density.

          if (dim <= this->m_lat->getOrder()) {
            if (this->m_dualF) // dual basis
              this->m_normalizer->setLogDensity(Dbl( - dim * log(this->m_lat->getModulo()) ));
            else // primal basis
              this->m_normalizer->setLogDensity(Dbl( dim * log(this->m_lat->getModulo())  ));
          }
          //std::cout << "density APRES = " << exp(this->m_normalizer->getLogDensity()) << std::endl;
          //std::cout << "this->m_normalizer->getBound = " << this->m_normalizer->getBound(dim) << std::endl;

          if (this->m_lat->getNorm () == LatticeTester::L2NORM) {

            if (!this->m_dualF) { // in case we work with rescaled values
              double normalizer = this->m_normalizer->getBound(dim);
              this->m_merit[dim] = log(temp) - 2*log(normalizer) - 2*log(mr);
              this->m_merit[dim] = exp(this->m_merit[dim]);
            } else { // general case
              double normalizer = this->m_normalizer->getBound(dim);
              this->m_merit[dim] = log(temp) - 2*log(normalizer);
              this->m_merit[dim] = exp(this->m_merit[dim]);
            }

          } else if (this->m_lat->getNorm () == LatticeTester::L1NORM) {

            if (!this->m_dualF) { // in case we work with rescaled values
              double normalizer = this->m_normalizer->getBound(dim);
              this->m_merit[dim] = log(temp) - log(normalizer) - log(mr);
              this->m_merit[dim] = exp(this->m_merit[dim]);
            } else { // general case
              double normalizer = this->m_normalizer->getBound(dim);
              this->m_merit[dim] = log(temp) - log(normalizer);
              this->m_merit[dim] = exp(this->m_merit[dim]);
            }

          } else
            this->m_merit[dim] = temp;

          this->m_merit[dim] /= weight;


          //std::cout << "this->m_merit[" << dim << "] = " << this->m_merit[dim] << std::endl;

          // Si on sait deja que ce gen. ne pourra etre retenu,
          // on le rejette tout de suite et on arrete le test.
          if ((this->m_maxAllDimFlag && this->m_merit[dim] < minVal[toDim])
              || this->m_merit[dim] < minVal[dim]) {
            //    this->m_merit[dim] = 0.0;
            return false;
          }
        }
        else
          this->m_merit[dim] /= weight;

        prepAndDisp (dim);

        if (dim == toDim)
          break;

        if (this->m_dualF)
          this->m_lat->dualize ();
        this->m_lat->incDim ();
        red = LatticeTester::Reducer<Int, Int, Dbl, Dbl>(*this->m_lat);
      }

      return true;
    }

  //===========================================================================

  template<typename Int, typename Dbl>
    bool LatTestSpectral<Int, Dbl>::test (int fromDim, int toDim, double minVal[], const double* weights)
    {

      //BOOST_DISPLAY
      //boost::progress_display show_progress(toDim-fromDim+1);

      this->m_merit.setDim(toDim);
      this->m_fromDim = fromDim;
      this->m_toDim = toDim;
      init ();
      this->resetFromDim (this->m_lat->getOrder (), fromDim);
      if (this->m_normalizer->getNorm () != this->m_lat->getNorm ()) {
        std::cout << "SpectralTest: conflict between NormType and Normalizer" <<
          std::endl;
        exit (EXIT_FAILURE);
      }

      double mr;
      NTL::conv (mr, this->m_lat->getModulo ());
      double temp;

      //Dbl te;
      //double lgvv1 = 0.0;
      // "unused variables" according to the compiler

      while (this->m_lat->getDim () < fromDim)
        this->m_lat->incDim ();

      LatticeTester::Reducer<Int, Int, Dbl, Dbl> red (*this->m_lat);

      if (this->m_S2toL2[fromDim] <= 0.0)
        initLowerBoundL2 (fromDim, toDim);
      setLowerBoundL2 (minVal[toDim], weights);   // same S2 for all dim
      red.setBoundL2 (this->m_boundL2, fromDim, toDim);

      while (true) {

        //BOOST_DISPLAY
        //++show_progress;

        if (this->m_dualF)
          this->m_lat->dualize ();

        int dim = this->m_lat->getDim ();

        // pre-reduction step before BB with default parameters
        red.redBKZ(0.999999, 10, LatticeTester::QUADRUPLE, dim);
        //PW_TODO: hard coding?

        //std::cout << "dual basis = \n" << this->m_lat->getBasis() << std::endl;
        //std::cout << "primal basis = \n" << this->m_lat->getDualBasis() << std::endl;

        if (red.shortestVector (this->m_lat->getNorm ())) {

          // Calcul de D2. Pour Norm # L2NORM, suppose que VV est a jour.
          if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
            this->m_lat->updateScalL2Norm (0);
            NTL::conv (temp, this->m_lat->getVecNorm (0));

            // to work with (0,1) variables otherwise we get the rescaled values
            /*if (!this->m_dualF) {
              NTL::conv(te, temp);
              Dbl m2;
              NTL::conv(m2, this->m_lat->getModulo ());
              m2 = m2*m2;
              te = te / m2;
              NTL::conv(temp, te);
              }*/

          } else {
            NTL::conv (temp, red.getMinLength ());

            // to work with (0,1) variables otherwise we get the rescaled values
            /*if (!this->m_dualF)
              temp = temp / mr;
              */
          }
          if (3 == this->m_detailF) {
            this->dispatchLatUpdate(*this->m_lat);
            /*if (this->m_dualF) {
              dispatchLatUpdate (this->m_lat->getDualBasis());
              dispatchLatUpdate (this->m_lat->getBasis());
              } else {
              dispatchLatUpdate (this->m_lat->getBasis());
              dispatchLatUpdate (this->m_lat->getDualBasis());
              }*/
          }
          // weight factor:
          // require a higher merit for more important projections by dividing
          // their merit by the weight of the projection
          double weight = weights ? weights[dim] : 1.0;

          this->m_merit.getMerit (dim) = temp / weight;

          if (dim <= toDim) { // Calcul de S2.

            //std::cout << "dim = " << dim << std::endl;
            //std::cout << "order = " << this->m_lat->getOrder() << std::endl;
            //std::cout << "density AVANT = " << exp(this->m_normalizer->getLogDensity()) << std::endl;

            //updating value of matrix density.

            if (dim <= this->m_lat->getOrder()) {
              if (this->m_dualF) // dual basis
                this->m_normalizer->setLogDensity(Dbl( - dim * log(this->m_lat->getModulo()) ));
              else // primal basis
                this->m_normalizer->setLogDensity(Dbl( dim * log(this->m_lat->getModulo())  ));
            }

            //std::cout << "density APRES = " << exp(this->m_normalizer->getLogDensity()) << std::endl;
            //std::cout << "this->m_normalizer->getBound = " << this->m_normalizer->getBound(dim) << std::endl;

            if (this->m_lat->getNorm () == LatticeTester::L2NORM) {

              if (!this->m_dualF) { // in case we work with rescaled values
                double normalizer = this->m_normalizer->getBound(dim);
                this->m_merit[dim] = log(temp) - 2*log(normalizer) - 2*log(mr);
                this->m_merit[dim] = exp(this->m_merit[dim]);
              } else { // general case
                double normalizer = this->m_normalizer->getBound(dim);
                this->m_merit[dim] = log(temp) - 2*log(normalizer);
                this->m_merit[dim] = exp(this->m_merit[dim]);
              }

            } else if (this->m_lat->getNorm () == LatticeTester::L1NORM) {

              if (!this->m_dualF) { // in case we work with rescaled values
                double normalizer = this->m_normalizer->getBound(dim);
                this->m_merit[dim] = log(temp) - log(normalizer) - log(mr);
                this->m_merit[dim] = exp(this->m_merit[dim]);
              } else { // general case
                double normalizer = this->m_normalizer->getBound(dim);
                this->m_merit[dim] = log(temp) - log(normalizer);
                this->m_merit[dim] = exp(this->m_merit[dim]);
              }

            } else
              this->m_merit[dim] = temp;

            this->m_merit[dim] /= weight;


            std::cout << "this->m_merit[" << dim << "] = " << this->m_merit[dim] << std::endl;

            // Si on sait deja que ce gen. ne pourra etre retenu,
            // on le rejette tout de suite et on arrete le test.
            if ((this->m_maxAllDimFlag && this->m_merit[dim] < minVal[toDim])
                || this->m_merit[dim] < minVal[dim]) {
              //    this->m_merit[dim] = 0.0;
              return false;
            }
          }
          else
            this->m_merit[dim] /= weight;

          prepAndDisp (dim);

        } else {
          this->m_merit[dim] = -1.0;
          return false;
        }

        if (dim == toDim)
          break;

        if (this->m_dualF)
          this->m_lat->dualize ();
        this->m_lat->incDim ();
        red = LatticeTester::Reducer<Int, Int, Dbl, Dbl>(*this->m_lat);
      }

      return true;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void LatTestSpectral<Int, Dbl>::init ()
    {
      this->m_merit.setDim(this->m_toDim);
      this->m_boundL2.SetLength (this->m_toDim+1);
      this->m_S2toL2 = new double[this->m_toDim+1];
      LatticeTester::SetZero (this->m_S2toL2, this->m_toDim+1);
      const int N = 3;
      std::string header[N];
      if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
        if (this->m_invertF)
          header[0] = "d_t";
        else
          header[0] = "l_t";
        header[1] = "S_t";
        header[2] = "Cumul CPU t(sec)";
        this->dispatchTestInit ("SPECTRAL", header, N);

      } else {
        if (this->m_invertF)
          header[0] = "d_t";
        else
          header[0] = "N_t";
        //      header[1] = "N_t^*";
        header[1] = "S_t";
        header[2] = "Cumul CPU t(sec)";
        this->dispatchTestInit ("SPECTRAL", header, N);
      }


      this->timer.init ();
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void LatTestSpectral<Int, Dbl>::prepAndDisp (int dim)
    {
      const int N = 3;
      double results[N];

      if (2 == this->m_detailF)
        this->dispatchLatUpdate (*this->m_lat);
      else if (1 == this->m_detailF)
        this->dispatchLatUpdate (*this->m_lat, 0);

      if (this->m_lat->getNorm () == LatticeTester::L2NORM) {
        results[0] = sqrt (this->m_merit.getMerit (dim));    // L_t
        if (this->m_invertF)
          results[0] = 1.0 / results[0];   // d_t = 1/L_t
        results[1] = sqrt (this->m_merit[dim]);

        results[2] = this->timer.val (Chrono::SEC);
        this->dispatchResultUpdate (results, N);

      } else {   // LatticeTester::L1NORM
        results[0] = this->m_merit.getMerit (dim);
        if (this->m_invertF)
          results[0] = 1.0 / results[0];   // d_t = 1/N_t

#if 0
        // Calcule la valeur de N_t^*
        double x;
        NTL::conv (x, this->m_lat->getM());
        x = pow (x, (double) this->m_lat->getOrder());
        double y;
        if (this->m_dualF) {
          y = floor (pow (Factorial (dim) * x, 1.0 / dim));
        } else {
          x = 1.0/x;
          y = pow (Factorial (dim) * x, 1.0 / dim);
        }
        results[1] = y;
#endif

        results[1] = this->m_merit[dim];
        results[2] = this->timer.val (Chrono::SEC);
        this->dispatchResultUpdate (results, N);
      }
    }

  extern template class LatTestSpectral<std::int64_t, double>;
  extern template class LatTestSpectral<NTL::ZZ, double>;
  extern template class LatTestSpectral<NTL::ZZ, NTL::RR>;

}
#endif

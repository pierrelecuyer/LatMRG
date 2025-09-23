#ifndef LATMRG_MIXMAXMMRG_H
#define LATMRG_MIXMAXMMRG_H

#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "latticetester/Util.h"

namespace LatMRG {

  /**
   * This class is used to manipulate easily MMRG of Mixmax types as
   * described by Savvidy [CITE PAPER].
   */
  template<typename Int>
    class MixmaxMMRG {
      private:
        typedef NTL::matrix<Int> IntMat;

      public:
        /**
         * Constructor for the four-parameters family.
         */
        MixmaxMMRG(Int modulus, int N, Int & s, Int & m, Int & b);

        /**
         * Constructor for the three-parameters family.
         * The four-parameter family reduces to the three-parameter family with
         * \f$b=2-2m\f$.
         */
        MixmaxMMRG(Int modulus, int N, Int & s, Int & m);

        /**
         * Constructor for the two-parameters family.
         * The three-parameter family reduces to the two-parameter family with
         * \f$m=1\f$.
         */
        MixmaxMMRG(Int modulus, int N, Int & s);

        /**
         * Destructor.
         */
        ~MixmaxMMRG();

        /**
         * Returns the order of the MMRG matrix.
         */
        int getOrder() { return m_order; }

        /**
         * Returns the matrix of the MMRG
         */
        IntMat getMatrix() { return m_A; }

      private:
        /**
         * This method build the matrix for the MMRG recurrence according
         * to the selected parameters for the four-parameters family of generators.
         */
        void buildMatrix(Int modulus, int N, Int & s, Int & m, Int & b);

        /**
         * The modulus used for the MMRG.
         */
        Int m_modulus;

        /**
         * The order of the MMRG.
         */
        int m_order;

        /**
         * The 'magic number' parameter of the MMRG matrix.
         */
        Int m_parameter1;

        /**
         * Second parameter of the MMRG matrix.
         */
        Int m_parameter2;

        /**
         * Third parameter of the MMRG matrix.
         */
        Int m_parameter3;

        /**
         * The matrix used for the MMRG recurrence.
         */
        IntMat m_A;

    };

  //===========================================================================

  template<typename Int>
  MixmaxMMRG<Int>::MixmaxMMRG(Int modulus, int N, Int & s, Int & m, Int & b)
  {
    m_modulus = modulus;
    m_order = N;
    m_parameter1 = s;
    m_parameter2 = m;
    m_parameter3 = b;
    buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

  }

  //===========================================================================

  template<typename Int>
  MixmaxMMRG<Int>::MixmaxMMRG(Int modulus, int N, Int & s, Int & m)
  {
    m_modulus = modulus;
    m_order = N;
    m_parameter1 = s;
    m_parameter2 = m;
    m_parameter3 = 2 - 2*m;
    buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

  }

  //===========================================================================

  template<typename Int>
  MixmaxMMRG<Int>::MixmaxMMRG(Int modulus, int N, Int & s)
  {
    m_modulus = modulus;
    m_order = N;
    m_parameter1 = s;
    m_parameter2 = 1;
    m_parameter3 = 0; // =2-2*m
    buildMatrix(m_modulus, m_order, m_parameter1, m_parameter2, m_parameter3);

  }

  //===========================================================================

  template<typename Int>
  MixmaxMMRG<Int>::~MixmaxMMRG()
  {}

  //===========================================================================

  template<typename Int>
  void MixmaxMMRG<Int>::buildMatrix(Int modulus, int N, Int & s, Int & m, Int & b)
  {
    m_A.kill();
    m_A.resize(N,N);
    for (int j = 1; j < N; j ++) {
      for (int i = j+1; i < N; i++) {
        m_A[i][j] = (i-j+2) * m + b;
        LatticeTester::Modulo(m_A[i][j], modulus, m_A[i][j]);
      }
    }
    for (int i = 0; i < N; i ++)
      m_A[i][0] = 1;
    for (int i = 1; i < N; i++)
      m_A[i][i] = 2;
    for (int i = 0; i < N; i ++) {
      for (int j = i+1; j < N; j ++)
        m_A[i][j] = 1;
    }
    m_A[2][1] += s;
    LatticeTester::Modulo(m_A[2][1], modulus, m_A[2][1]);
  }

  template class MixmaxMMRG<std::int64_t>;
  template class MixmaxMMRG<NTL::ZZ>;

} // end namespace LatMRG

#endif

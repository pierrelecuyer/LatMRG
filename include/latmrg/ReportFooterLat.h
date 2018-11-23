#ifndef REPORTFOOTERLAT_H
#define REPORTFOOTERLAT_H

#include "latticetester/Writer.h"
#include "latticetester/Const.h"

#include "latmrg/ReportFooter.h"
#include "latmrg/LatticeTest.h"
#include "latmrg/Merit.h"


namespace LatMRG {

  /**
   * This class is an implementation of the `ReportFooter` abstract class for
   * the program `lat*`.
   *
   */
  template<typename Int, typename Dbl>
    class ReportFooterLat : public ReportFooter<Int> {
      public:

        /**
         * Constructor. `writer` is the writing engine used to write the report
         * footer. `test` is the lattice test thas was performed and for which the
         * results are to be written.
         */
        ReportFooterLat (LatticeTester::Writer<Int> * writer, LatticeTest<Int, Dbl> * test = 0): ReportFooter<Int>(writer) {
          this->m_test = test;
        }


        /**
         * `test` is the lattice test thas was performed and for which the
         * results are to be written.
         */
        void setLatticeTest (LatticeTest<Int, Dbl> * test) { this->m_test = test; }

        /**
         * Defined in abstract class `ReportFooter`.
         */
        void printFooter()
        {
          int dimMin = this->m_test->getMinDim ();
          int dimMax = this->m_test->getMaxDim ();
          int dimWorst;
          // Get minimal merit value and dimension where it occurs
          double S = this->m_test->getMerit().getST (dimMin, dimMax, dimWorst);
          if (this->m_test->getLattice()->getNorm() == LatticeTester::L2NORM)
            S = sqrt(S);
          this->m_writer->newLine ();
          this->m_writer->writeString (" Min merit:   S_");
          this->m_writer->writeInt (dimWorst);
          this->m_writer->writeString (" = ");
          this->m_writer->writeDouble (S);
          this->m_writer->newLine ();
          this->m_writer->newLine ();
        }
      private:

        /**
         * Pointer to the final lattice test which was performed.
         */
        LatticeTest<Int, Dbl> * m_test;
    };

}
#endif

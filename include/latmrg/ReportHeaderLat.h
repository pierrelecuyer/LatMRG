#ifndef	REPORTHEADERLAT_H
#define	REPORTHEADERLAT_H

#include <string>

#include "latticetester/Const.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Writer.h"

#include "latmrg/ReportHeader.h"
#include "latmrg/LatConfig.h"

namespace LatMRG {

  /**
   * This class is an implementation of the `ReportHeader` abstract class for
   * the programs `lat*`. It prints the configuration of the test launched, the
   * MRG or MWC components, and the combined MRG if applicable.
   *
   */
  template<typename Int, typename Dbl>
    class ReportHeaderLat : public ReportHeader<Int> {
      public:

        /**
         * Constructor. `writer` is the writing engine used to write the report
         * header. `config` is the configuration of the lattice test to be performed,
         * which is populated from an instance of `ParamReaderLat`. `lattice` is the
         * final MRG lattice on which the lattice test will be performed. It can be
         * the result of a combination of MRG components, a MWC transformed into a
         * MRG, and so on.
         */
        ReportHeaderLat (LatticeTester::Writer<Int> *writer, LatConfig<Int> *config,
            LatticeTester::IntLattice<Int, Int, Dbl, Dbl> *lattice);

        /**
         * Does the actual writing of the report header with the `Writer`
         * passed to the constructor.
         */
        void printHeader();
      private:

        /**
         * Pointer to the configuration of the lattice test.
         */
        LatConfig<Int>* m_config;

        /**
         * Pointer to the final MRG lattice on which the lattice test will be
         * performed.
         */
        LatticeTester::IntLattice<Int, Int, Dbl, Dbl>* m_lattice;
    };

  //============================================================================

  template<typename Int, typename Dbl>
    ReportHeaderLat<Int, Dbl>::ReportHeaderLat (LatticeTester::Writer<Int> * writer, LatConfig<Int> * config,
        LatticeTester::IntLattice<Int, Int, Dbl, Dbl> * lattice): ReportHeader<Int> (writer)
  {
    m_config = config;
    m_lattice = lattice;
  }


  template<typename Int, typename Dbl>
    void ReportHeaderLat<Int, Dbl>::printHeader ()
    {
      this->m_writer->newLine ();
      this->m_writer->writeString (
          "===========================================================================");
      this->m_writer->newLine ();
      this->m_writer->writeString ("Generator's components:");
      this->m_writer->newLine ();
      this->m_writer->writeString ("Number of components: ");
      this->m_writer->writeInt (m_config->J);
      this->m_writer->newLine ();

      // Écriture des informations sur les composantes

      this->m_writer->beginTabbedSection ();
      this->m_writer->addTab ();

      for (int i = 0; i < m_config->J; i++) {
        if (m_config->J > 1) {
          this->m_writer->newLine ();
          this->m_writer->writeString ("Component ");
          this->m_writer->writeInt (i + 1);
          this->m_writer->newLine ();

        }

        this->m_writer->writeString ("   Lattice type: ");
        this->m_writer->writeString (toStringGen (m_config->genType[i]));
        this->m_writer->newLine ();

        this->m_writer->writeString ("   m = ");
        this->m_writer->writeIntScal (m_config->comp[i]->getM ());

        this->m_writer->newLine ();

        if (m_config->genType[i] == MMRG) {
          this->m_writer->writeString ("   order = ");
          this->m_writer->writeInt (m_config->comp[i]->k);
          this->m_writer->newLine ();

          this->m_writer->writeString ("   generator matrix = ");
          this->m_writer->newLine ();
          this->m_writer->writeMMat(m_config->comp[i]->A);


        } else { 
          this->m_writer->writeString ("   k = ");
          this->m_writer->writeInt (m_config->comp[i]->k);
          this->m_writer->newLine ();

          if (m_config->comp[i]->k < 13) {
            // Write all coefficients of a
            for (int ii = 0; ii < m_config->comp[i]->k; ii++) {
              this->m_writer->writeString ("      a_");
              this->m_writer->writeInt (ii);
              this->m_writer->writeString (" = ");
              this->m_writer->writeIntScal (m_config->comp[i]->a[ii]);
              this->m_writer->newLine ();
            }
          } else {
            // Write only the non-zero coefficients of a
            for (int ii = 0; ii < m_config->comp[i]->k; ii++) {
              if (0 != m_config->comp[i]->a[ii]) {
                this->m_writer->writeString ("      a_");
                this->m_writer->writeInt (ii);
                this->m_writer->writeString (" = ");
                this->m_writer->writeIntScal (m_config->comp[i]->a[ii]);
                this->m_writer->newLine ();
              }
            }
            this->m_writer->writeString ("         All other a_j are  0");
            this->m_writer->newLine ();
          }
        }

        if (m_config->J > 1) {
          this->m_writer->writeString ("   nj = ");
          this->m_writer->writeIntScal (m_config->comp[i]->nj);
          this->m_writer->newLine ();
          this->m_writer->writeString ("   rho = ");
          this->m_writer->writeIntScal (m_config->comp[i]->rho);
          this->m_writer->newLine ();
        }
      }

      if (m_config->J > 1) {
        this->m_writer->newLine ();
        this->m_writer->writeString ("Combined MRG:");
        this->m_writer->newLine ();
        this->m_writer->writeString ("   m = ");
        this->m_writer->writeIntScal (m_lattice->getModulo());
        this->m_writer->newLine ();
        this->m_writer->writeString ("   k = ");
        this->m_writer->writeInt (m_lattice->getOrder ());
        this->m_writer->newLine ();

        for (int i = 1; i <= m_lattice->getOrder (); i++) {
          this->m_writer->writeString ("      a_");
          this->m_writer->writeInt (i);
          this->m_writer->writeString (" = ");
          //this->m_writer->writeIntScal (m_lattice->getCoef ()[i]);
          this->m_writer->newLine ();
        }
        this->m_writer->writeString ("   rho = ");
        //this->m_writer->writeIntScal (m_lattice->getRho ());
        this->m_writer->newLine ();
        this->m_writer->writeString ("   lossRho = ");
        //this->m_writer->writeIntScal (m_lattice->getLossRho ());
        this->m_writer->newLine ();
      }

      this->m_writer->newLine ();
      this->m_writer->endTabbedSection ();

      this->m_writer->writeString ("Test: ");
      this->m_writer->writeString (toStringCriterion (m_config->criter));
      this->m_writer->newLine ();
      if (m_config->criter == LatticeTester::PALPHA) {
        this->m_writer->writeString ("CalcType: ");
        this->m_writer->writeString (toStringCalc (m_config->calcPalpha));
        this->m_writer->newLine ();
      }

      if (m_config->criter == LatticeTester::SPECTRAL) {
        this->m_writer->writeString ("Norm: ");
        this->m_writer->writeString (toStringNorm (m_config->norm));
        this->m_writer->newLine ();
        this->m_writer->writeString ("Normalisation: ");
        this->m_writer->writeString (toStringNorma (m_config->norma));
        this->m_writer->newLine ();
      }

      if (m_config->criter == LatticeTester::PALPHA) {
        this->m_writer->writeString ("Prime m: ");
        this->m_writer->writeBool (m_config->primeM);
        this->m_writer->newLine ();
        this->m_writer->writeString ("Max period: ");
        this->m_writer->writeBool (m_config->maxPeriod);
        this->m_writer->newLine ();
        this->m_writer->writeString ("Alpha: ");
        this->m_writer->writeDouble (m_config->alpha);
        this->m_writer->newLine ();
        this->m_writer->writeString ("Seed: ");
        this->m_writer->writeInt (m_config->seed);
        this->m_writer->newLine ();
        this->m_writer->writeString ("Beta: { ");
        this->m_writer->writeDouble (m_config->Beta[0]);
        for (int i = 0; i <= m_config->td[1]; i++) {
          this->m_writer->writeString (", ");
          this->m_writer->writeDouble (m_config->Beta[i]);
        }
        this->m_writer->writeString (" }");

      } else {
        this->m_writer->writeString ("Lattice Type: ");
        this->m_writer->writeString (toStringLattice (m_config->latType));
        this->m_writer->newLine ();
        if (m_config->dualF)
          this->m_writer->writeString ("Lattice: DUAL");
        else
          this->m_writer->writeString ("Lattice: PRIMAL");
        this->m_writer->newLine ();
        this->m_writer->writeString ("Dimensions td:");
        for (int i = 0; i <= m_config->d; i++) {
          this->m_writer->writeString (" ");
          this->m_writer->writeInt (m_config->td[i]);
        }
        this->m_writer->newLine ();
      }

      if (m_config->lacunary) {
        int dim = m_config->td[1];
        this->m_writer->writeString ("dim = ");
        this->m_writer->writeInt (dim);
        if (m_config->genType[0] == MMRG) {

          if (m_config->lacunaryType == SUBVECTOR)
            this->m_writer->writeString ("\n\nLacunary indices for each vector = {   ");
          else if (m_config->lacunaryType == ARBITRARYINDICES)
            this->m_writer->writeString ("\n\nLacunary indices = {   ");

          Int pre;
          NTL::conv (pre, -9);
          for (int i = 0; i < m_config->numberLacIndices; i++) {
            Int r = m_config->Lac[i];
            if (pre < r - 1)
              this->m_writer->writeString ("\n   ");
            this->m_writer->writeIntScal (r);
            this->m_writer->writeString ("    ");
            pre = r;
          }

        } else {
          this->m_writer->writeString ("\n\nLacunary = {   ");
          // Print indices by group of s
          Int pre;
          NTL::conv (pre, -9);
          for (int i = 0; i < dim; i++) {
            Int r = m_config->Lac[i];
            if (pre < r - 1)
              this->m_writer->writeString ("\n   ");
            this->m_writer->writeIntScal (r);
            this->m_writer->writeString ("    ");
            pre = r;
          }
        }
        this->m_writer->writeString ("\n}\n");
      }

      this->m_writer->newLine ();
      this->m_writer->newLine ();
    }

}
#endif

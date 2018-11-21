#ifndef REPORTLAT_H
#define REPORTLAT_H

#include "latticetester/Writer.h"
#include "latticetester/Types.h"

#include "latmrg/LatticeTestObserver.h"
#include "latmrg/TableColumnImpl.h"
#include "latmrg/ReportHeader.h"
#include "latmrg/ReportFooter.h"
#include "latmrg/DoubleFormatter.h"
#include "latmrg/FormatterImpl.h"
#include "latmrg/LatConfig.h"
#include "latmrg/Table.h"

#include <string>
#include <iostream>


namespace LatMRG {

  /**
   * This class formats and prints the actual report for a lattice test for the
   * program `lat*`. It implements the interface `LatticeTestObserver` to be
   * able to receive information and results from the lattice test.
   */
  template<typename Int, typename Dbl>
    class ReportLat : public LatticeTestObserver<Int, Dbl> {
      public:

        /**
         * Constructor.
         */
        ReportLat (LatticeTester::Writer<Int>* writer, LatConfig<Int>* config,
            ReportHeader<Int>* header, ReportFooter<Int>* footer);

        /**
         * Destructor.
         */
        ~ReportLat();

        /**
         * Prints the header using the `ReportHeader` passed to the
         * constructor.
         */
        void printHeader();

        /**
         * Prints the footer using the `ReportFooter` passed to the
         * constructor.
         */
        void printFooter();

        /**
         * Prints the table of results obtained from the successive calls of
         * `resultUpdate`. If more than one tests are performed, the results
         * will be concatenated in the same table.
         */
        void printTable();

        /**
         * Defined in interface `LatticeTestObserver`. Prints the base directly
         * in the report.
         */
        void latUpdate (LatticeTester::IntLattice<Int, Int, Dbl, Dbl> & lat);

        /**
         * Defined in interface `LatticeTestObserver`. Prints basis vector
         * \f$V[i]\f$ directly in the report.
         */
        void latUpdate (LatticeTester::IntLattice<Int, Int, Dbl, Dbl> & lat, int i);

        /**
         * Defined in interface `LatticeTestObserver`. The results are stacked
         * in the internal table and will be printed upon a call to
         * `printTable`.
         */
        void resultUpdate (double[], int);

        /**
         * Defined in interface `LatticeTestObserver`. The columns in the
         * internal table are set up to be able to receive results from calls
         * to `resultUpdate`.
         */
        void testInit (const std::string &, std::string[], int);

        /**
         * Defined in interface `LatticeTestObserver`. Indicates a successful
         * test.
         */
        void testCompleted();

        /**
         * Defined in interface `LatticeTestObserver`. Indicates a failed test.
         * An error message is printed in the report.
         */
        void testFailed (int);

        /**
         * Returns the writing engine used in this class
         */
        LatticeTester::Writer<Int> * getWriter() { return m_writer; }
      private:

        /**
         * Writing engine used to print the report.
         */
        LatticeTester::Writer<Int> * m_writer;

        /**
         * Report header used to print the report.
         */
        ReportHeader<Int> * m_header;

        /**
         * Report footer used to print the report.
         */
        ReportFooter<Int> * m_footer;

        /**
         * Indicates the index of the first column in the results table to
         * insert results for the current test.
         */
        int m_base_col;

        /**
         * Indicates the index of the first column in the results table to
         * insert results for the next test.
         */
        int m_next_base;

        /**
         * The internal table that contains the results of lattice test(s).
         */
        Table m_results;

        /**
         * The configuration of the lattice test(s).
         */
        LatConfig<Int> * m_config;

        /**
         * Formatter used to format the results in the table when writing the
         * report.
         */
        DoubleFormatter m_dFormat;

        /**
         * Formatter used to format the dimensions (first column) in the table
         * when writing the report.
         */
        FormatterImpl<int> m_iFormat;
    };

  template<typename Int, typename Dbl>
    ReportLat<Int, Dbl>::ReportLat(LatticeTester::Writer<Int>* writer, LatConfig<Int>* config,
        ReportHeader<Int>* header, ReportFooter<Int>* footer): m_dFormat(6)
  {
    m_config = config;
    m_writer = writer;
    m_header = header;
    m_footer = footer;

    TableColumnImpl<int>* dims = new TableColumnImpl<int>(&m_iFormat, "t");

    m_results.add(dims);

    for (int i = m_config->td[0]; i <= m_config->td[1]; i++) {
      dims->addValue(static_cast<void*>(&i));
    }

    m_base_col = 1;
  }

  template<typename Int, typename Dbl>
    ReportLat<Int, Dbl>::~ReportLat()
    {
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::printHeader()
    {
      m_header->printHeader();
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::printFooter()
    {
      m_footer->printFooter();
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::printTable()
    {
      std::ostream* stream(&this->getWriter()->getStream());
      Table data = m_results;
      const std::string pos = "llllllllll";
      int margins = 5;
      int t = data.size ();       // Nb de colonne
      int n = data.getHeight () + 1; // Nb de ligne

      int max_length[t];
      int nb_lig_cell[t][n];
      int nb_lig[n];
      std::string::size_type curpos;
      int curpos_b;
      int lig;
      int total_length = 0;

      memset (max_length, 0, t * sizeof (int));
      memset (nb_lig, 0, n * sizeof (int));

      std::string temp;

      // On calcule le nombre de sous ligne dans chaque case et la longueur
      // maximale de chaque colonne.

      for (int i = 0; i < t; i++) {
        for (int j = -1; j < data[i]->size (); j++) {
          lig = 0;
          curpos = 0;
          curpos_b = -1;

          temp = data[i]->getFormattedValue (j);

          while ((curpos = temp.find ("\n", curpos)) != std::string::npos) {
            int tmp = (int) curpos - curpos_b - 1;
            max_length[i] = std::max (max_length[i], tmp);
            lig++;
            curpos_b = (int) curpos;
            curpos++;
          }
          if (curpos_b == -1) {
            curpos_b = 0;
          }

          max_length[i] = std::max (max_length[i],
              static_cast <int>(temp.length () - curpos_b));
          nb_lig[j + 1] = std::max (nb_lig[j + 1], ++lig);
          nb_lig_cell[i][j + 1] = lig;
        }

        total_length += max_length[i] + 2 * margins;
      }

      int nb_sub_lig = 0;

      for (int i = 0; i < n; i++) {
        nb_sub_lig += nb_lig[i];
      }

      int k, offset;



      //string sub_data[t][nb_sub_lig];
      // this is not compiling with old versions of clang

      std::string **sub_data = new std::string* [t];
      for (int i = 0; i < t; ++i)
        sub_data[i] = new std::string [nb_sub_lig];
      // better could be : string *sub_data = new string [t * nb_sub_lig]



      // On décompose les chaines en sous chaines

      for (int i = 0; i < t; i++) {
        offset = 0;
        for (int j = -1; j < data[i]->size (); j++) {
          curpos = 0;
          curpos_b = 0;

          // padding with empty std::strings
          for (k = 0; k < nb_lig[j + 1] - nb_lig_cell[i][j + 1]; k++) {
            sub_data[i][j + offset + k + 1] = "";
          }

          temp = data[i]->getFormattedValue (j);

          while ((curpos = temp.find ("\n", curpos)) != std::string::npos) {
            sub_data[i][j + offset + k + 1] =
              temp.substr (curpos_b, curpos);
            k++;
            curpos_b = (int) curpos;
            curpos++;
          }
          if (curpos_b > 0) {
            curpos_b++;
          }
          sub_data[i][j + offset + k + 1] =
            temp.substr (curpos_b, temp.length () - curpos_b);
          sub_data[i][j + offset + k + 1] = temp.substr (curpos_b, curpos);
          offset += k;
        }
      }

      // On dessine le tableau, finalement ...

      int length;
      char align;
      for (int i = 0; i < nb_sub_lig; i++) {
        for (int j = 0; j < t; j++) {

          align = pos[j];
          length = (int) sub_data[j][i].length ();

          switch (align) {
            case 'c':
              offset = (max_length[j] >> 1) - (length >> 1);
              break;
            case 'r':
              offset = max_length[j] - length;
              break;
            default:               // case l
              offset = 0;

          }

          for (int k = 0; k < margins + offset; k++) {
            *stream << " ";
          }

          *stream << sub_data[j][i];

          for (int k = 0; k < max_length[j] - length - offset + margins; k++) {
            *stream << " ";
          }
        }
        if (i == nb_lig[0] - 1) {
          *stream << std::endl;
          for (int k = 0; k < total_length; k++) {
            *stream << "-";
          }
        }
        *stream << std::endl;
      }

      // deleting pointers
      for (int i = 0; i < t; ++i)
        delete [] sub_data[i];
      delete [] sub_data;

    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::latUpdate(LatticeTester::IntLattice<Int, Int, Dbl, Dbl> & lat)
    {
      m_writer->writeString (lat.toStringBasis());
      m_writer->writeString (lat.toStringDualBasis());
      //FIX ME
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::latUpdate(LatticeTester::IntLattice<Int, Int, Dbl, Dbl> & lat, int i)
    {
      m_writer->writeString (lat.toStringBasis());
      //FIX ME
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::resultUpdate(double results[], int n)
    {
      for (int i = 0; i < n; i++) {
        m_results[m_base_col + i]->addValue(static_cast<void*>(&results[i]));
      }
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::testInit(const std::string & test, std::string headers[], int n)
    {
      TableColumnImpl<double>* col;

      for (int i = 0; i < n; i++) {
        col = new TableColumnImpl<double>(&m_dFormat, headers[i]);
        m_results.add(col);
      }

      m_next_base = m_base_col + n;

      //   m_writer->writeString("Test: ");
      //   m_writer->writeString(test);
      //   m_writer->newLine();
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::testCompleted()
    {
      m_writer->newLine();
      m_writer->writeString("Test completed successfully");
      m_writer->newLine();

      m_base_col = m_next_base;
    }

  template<typename Int, typename Dbl>
    void ReportLat<Int, Dbl>::testFailed(int dim)
    {
      m_writer->newLine();
      m_writer->writeString("Test failed at dim: ");
      m_writer->writeInt(dim);
      m_writer->newLine();

      m_base_col = m_next_base;
    }

}
#endif

/*
   ReportLat.cc for ISO C++
   version 1.00

authors: Hicham Wahbi
Frédérik Rozon
Richard Simard
*/

#include "latmrg/ReportLat.h"
#include "latmrg/TableColumnImpl.h"

#include "latticetester/Writer.h"
#include "latticetester/Types.h"

using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

  ReportLat::ReportLat(Writer<MScal>* writer, LatConfig<MScal>* config, ReportHeader* header,
      ReportFooter* footer): m_dFormat(6)
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

  ReportLat::~ReportLat()
  {
  }

  void ReportLat::printHeader()
  {
    m_header->printHeader();
  }

  void ReportLat::printFooter()
  {
    m_footer->printFooter();
  }

  void ReportLat::printTable()
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
    string::size_type curpos;
    int curpos_b;
    int lig;
    int total_length = 0;

    memset (max_length, 0, t * sizeof (int));
    memset (nb_lig, 0, n * sizeof (int));

    string temp;

    // On calcule le nombre de sous ligne dans chaque case et la longueur
    // maximale de chaque colonne.

    for (int i = 0; i < t; i++) {
      for (int j = -1; j < data[i]->size (); j++) {
        lig = 0;
        curpos = 0;
        curpos_b = -1;

        temp = data[i]->getFormattedValue (j);

        while ((curpos = temp.find ("\n", curpos)) != string::npos) {
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

    string **sub_data = new string* [t];
    for (int i = 0; i < t; ++i)
      sub_data[i] = new string [nb_sub_lig];
    // better could be : string *sub_data = new string [t * nb_sub_lig]



    // On décompose les chaines en sous chaines

    for (int i = 0; i < t; i++) {
      offset = 0;
      for (int j = -1; j < data[i]->size (); j++) {
        curpos = 0;
        curpos_b = 0;

        // padding with empty strings
        for (k = 0; k < nb_lig[j + 1] - nb_lig_cell[i][j + 1]; k++) {
          sub_data[i][j + offset + k + 1] = "";
        }

        temp = data[i]->getFormattedValue (j);

        while ((curpos = temp.find ("\n", curpos)) != string::npos) {
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
        *stream << endl;
        for (int k = 0; k < total_length; k++) {
          *stream << "-";
        }
      }
      *stream << endl;
    }

    // deleting pointers
    for (int i = 0; i < t; ++i)
      delete [] sub_data[i];
    delete [] sub_data;

  }

  void ReportLat::latUpdate(IntLattice<MScal, BScal, NScal, RScal> & lat)
  {
    m_writer->writeString (lat.toStringBasis());
    m_writer->writeString (lat.toStringDualBasis());
    //FIX ME
  }

  void ReportLat::latUpdate(IntLattice<MScal, BScal, NScal, RScal> & lat, int i)
  {
    m_writer->writeString (lat.toStringBasis());
    //FIX ME
  }

  void ReportLat::resultUpdate(double results[], int n)
  {
    for (int i = 0; i < n; i++) {
      m_results[m_base_col + i]->addValue(static_cast<void*>(&results[i]));
    }
  }

  void ReportLat::testInit(const string & test, string headers[], int n)
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

  void ReportLat::testCompleted()
  {
    m_writer->newLine();
    m_writer->writeString("Test completed successfully");
    m_writer->newLine();

    m_base_col = m_next_base;
  }

  void ReportLat::testFailed(int dim)
  {
    m_writer->newLine();
    m_writer->writeString("Test failed at dim: ");
    m_writer->writeInt(dim);
    m_writer->newLine();

    m_base_col = m_next_base;
  }
}

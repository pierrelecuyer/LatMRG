#ifndef	REPORTHEADER_H
#define REPORTHEADER_H

#include "latticetester/Writer.h"

namespace LatMRG {

  /**
   * This is an abstract class that must implemented to print a header in
   * `ReportLat` or `ReportSeek`.
   *
   */
  template<typename Int>
  class ReportHeader {
    public:

      /**
       * Constructor. See the module `Writer` for more information.
       */
      ReportHeader (LatticeTester::Writer<Int> *writer) { m_writer = writer; }

      /**
       * Destructor.
       */
      virtual ~ReportHeader() {}

      /**
       * Writes the report header using the `Writer` passed to the
       * constructor.
       */
      virtual void printHeader() = 0;
    protected:

      /**
       * The `Writer` used to write the report header.
       */
      LatticeTester::Writer<Int>* m_writer;
  };

}
#endif

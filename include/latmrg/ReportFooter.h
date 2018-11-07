#ifndef	REPORTFOOTER_H
#define REPORTFOOTER_H

#include "latticetester/Writer.h"
#include "latticetester/Types.h"

namespace LatMRG {

  /**
   * This is an abstract class that must be implemented to print a footer in
   * `ReportLat` or `ReportSeek`.
   *
   */
  class ReportFooter {
    public:

      /**
       * Constructor. See the module `Writer` for more information.
       */
      ReportFooter (LatticeTester::Writer<MScal> *writer) {m_writer = writer;}

      /**
       * Destructor.
       */
      virtual ~ReportFooter() {}

      /**
       * Writes the report footer using the `Writer` passed to the
       * constructor.
       */
      virtual void printFooter() = 0;
    protected:

      /**
       * The `Writer` to be used to write the report footer.
       */
      LatticeTester::Writer<MScal>* m_writer;
  };

}
#endif

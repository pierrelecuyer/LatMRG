#ifndef DOUBLEFORMATTER_H
#define DOUBLEFORMATTER_H
#include "Formatter.h"
#include <string>


namespace LatMRG {

  /**
   * This class is an implementation of the `Formatter` interface to format
   * a `double` into a `string`.
   *
   */
  class DoubleFormatter : public Formatter {
    public:

      /**
       * Constructor initializing the object to format <tt>double</tt>’s with
       * `precision` decimals.
       */
      DoubleFormatter (int precision);

      /**
       * Formats the double `value` into a string. The argument `value` is
       * assumed to be a pointer to a `double`.
       */
      std::string format (void* value);
    private:

      /**
       * Precision used by the formatter.
       */
      int m_precision;
  };

}
#endif

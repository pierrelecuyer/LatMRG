#ifndef FORMATTERIMPL_H
#define FORMATTERIMPL_H
#include "Formatter.h"
#include <string>
#include <sstream>


namespace LatMRG {

  /**
   * This class is a template that implements the interface `Formatter`. It can
   * format any type with the `<<` operator defined.
   */

  template <typename T>
    class FormatterImpl: public Formatter {
      public:

        /**
         * Formats `value`, which is supposed to be of type `T`, into a string.
         */
        std::string format (void* value);
    };
  template <typename T>
    std::string FormatterImpl<T>::format(void* value) {
      T val;
      std::ostringstream outs;

      val = *static_cast<T*>(value);

      outs << val;

      return outs.str();

    }
}
#endif

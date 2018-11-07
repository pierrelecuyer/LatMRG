#ifndef FORMATTER_H
#define FORMATTER_H
#include <string>


namespace LatMRG {

  /**
   * This class is an interface that can be implemented to format different 
   * types into a `string`. It is used to format values in a `TableColumn` which
   * composes a `Table`. There is a general implementation of this class
   * available in the class `FormatterImpl`.
   *
   * \todo Since this class is implemented in a template class, is it necessary
   * to declare an interface?
   */
  class Formatter {
    public:

      /**
       * Destructor.
       */
      virtual ~Formatter() {}

      /**
       * Method that must implemented to format `value` into a string.
       */
      virtual std::string format (void *value) = 0;
  };

}
#endif

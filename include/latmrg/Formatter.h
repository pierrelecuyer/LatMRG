#ifndef FORMATTER_H
#define FORMATTER_H
#include <string>


namespace LatMRG {

  /**
   * This class is an interface that can be implemented to format diffenrent 
   * types into a `string`. It is used, therefore it must be implemented, to 
   * format values in a `TableColumn` which composes a `Table`.
   *
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


#ifndef LATMRGROUTINES_H
#define LATMRGROUTINES_H

#include "latmrg/LatConfig.h"

namespace LatMRG {

/**
 * This function allows computation of the shortest non-zero vector in a lattice,
 * according to the selected norm. Many parameters can bet set by the user, otherwise
 * the function work with default values.
 * Returns -1.0 if there was an error in Branch-and-Bound procedure. Return the length
 * of the shortest non-zero vector otherwise.
 */
std::vector<double> doLETest (LatConfig config);



} // end namespace LatticeTester

#endif
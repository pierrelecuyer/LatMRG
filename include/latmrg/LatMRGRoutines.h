
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
std::vector<double> ComputeFigureOfMerit (LatConfig& config);


void printResult(const std::vector<double> & result, const int & fromDim);

// une fonction qui calcule le shortest vector

// une fonction d'aide à la configuration de LatConfig
void initConfigSpectralTest(LatConfig& config);

void initConfigBeyerTest(LatConfig& config);

void initConfigPalphaTest(LatConfig& config);


} // end namespace LatticeTester

#endif
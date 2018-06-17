#ifndef MIXMAXMMRG_H
#define MIXMAXMMRG_H

#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "latticetester/Types.h"

namespace LatMRG {

/**
 * This class is used to manipulate easily MMRG of Mixmax types as
 * described by Savvidy [CITE PAPER].
 */

class MixmaxMMRG {

	public:
		/**
		 * Constructor for the four-parameters family.
		 */
		MixmaxMMRG(MScal modulus, int & N, MScal & s, MScal & m, MScal & b);

		/**
		 * Constructor for the three-parameters family. 
		 * The four-parameter family reduces to the three-parameter family with 
		 * \f$b=2-2m\f$.
		 */
		MixmaxMMRG(MScal modulus, int & N, MScal & s, MScal & m);

		/**
		 * Constructor for the two-parameters family.
		 * The three-parameter family reduces to the two-parameter family with 
		 * \f$m=1\f$.
		 */
		MixmaxMMRG(MScal modulus, int & N, MScal & s);

		/**
		 * Destructor.
		 */
		~MixmaxMMRG();

		/**
		 * Returns the order of the MMRG matrix.
		 */
		int getOrder() { return m_order; }

		/**
		 * Returns the matrix of the MMRG
		 */
		MMat getMatrix() { return m_A; }

	private:
		/**
		 * This method build the matrix for the MMRG recurrence according
		 * to the selected parameters for the four-parameters family of generators.
		 */
		void buildMatrix(MScal modulus, int & N, MScal & s, MScal & m, MScal & b);

		/**
		 * The modulus used for the MMRG.
		 */
		MScal m_modulus;

		/**
		 * The order of the MMRG.
		 */
		int m_order;

		/**
		 * The 'magic number' parameter of the MMRG matrix.
		 */
		 MScal m_parameter1;

		/**
		 * Second parameter of the MMRG matrix.
		 */
		 MScal m_parameter2;

		/**
		 * Third parameter of the MMRG matrix.
		 */
		 MScal m_parameter3;

		/**
		 * The matrix used for the MMRG recurrence.
		 */
		 MMat m_A;

};

} // end namespace LatticeTester

#endif

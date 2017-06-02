#include "latmrg/HyperplaneDistanceMerit.h"
#include "latticetester/Reducer.h"
#include "latmrg/Rank1Lattice.h"
#include "latticetester/NormaBestLat.h"

//!#define DEBUG
#undef DEBUG

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <cfloat>

#ifdef DEBUG
#include <iostream>
#endif


using namespace std;
using namespace LatticeTester;

namespace LatMRG
{

//===========================================================================

HyperplaneDistanceMerit::HyperplaneDistanceMerit (
   const Normalizer& normalizer, double power)
      : m_normalizer(normalizer), m_power(power)
{}


//===========================================================================

const std::string HyperplaneDistanceMerit::toString() const
{
   return "normalized maximum distance between successive hyperplanes";
}


//===========================================================================

double HyperplaneDistanceMerit::compute (const std::vector<long>& a, int n,
                        const std::vector<int>& projection) const
{
#ifdef DEBUG
   std::cout << "    computing spectral test in dimension "
             << projection.size() << std::endl;
   //! std::cout << " with generating vector [";
   //! for (std::vector<long>::const_iterator it = a.begin() + 1; it != a.end(); ++it)
   //!    std::cout << " " << *it;
   //! std::cout << " ]" << std::endl;
#endif

   if (m_normalizer.getNorm() != L2NORM)
      // this is the L2NORM implementation
      throw std::invalid_argument("norm of normalizer must be L2NORM");
   //if (projection.size() <= 1)
   //   throw std::invalid_argument("projection order must be >= 2");

   // projections of order 1 are the same for all lattices
   if (projection.size() <= 1)
      return 1.0;

   int projDim = projection.size(); // projection order

   // convert to NTL
   MVect big_a;
   big_a.SetLength (projDim);
   //big_a[0] = 0;
   for (int j = 0; j < projDim; j++)
      big_a[j] = a[projection[j]];

#ifdef DEBUG

   std::cout << "      projected vector: [";
   for (int j = 0; j < projDim; j++)
      std::cout << " " << big_a[j];
   std::cout << " ]" << std::endl;
#endif

   // convert to NTL
   MScal big_n;
   conv (big_n, n);

   // prepare lattice and basis reduction
   Rank1Lattice lattice (big_n, big_a, projDim, m_normalizer.getNorm());
   lattice.buildBasis (projDim);
   lattice.dualize ();
   Reducer reducer (lattice);

   if (!reducer.shortestVector(lattice.getNorm()))
      // reduction failed
      return DBL_MAX;

   // get length of shortest vector under L2NORM
   double sqlength;
   lattice.updateScalL2Norm(0);
   conv (sqlength, lattice.getVecNorm(0)); // square length

   double sqlength0 = m_normalizer.getGamma(projDim) * pow(n, 2.0 / projDim);
   // FIXME: use logarithms for big integers
   // FIXME: warn Richard about m_lgVolDual2 not being initialized in IntLattice through calcLgVolDual2()

   double merit = sqrt(sqlength0 / sqlength);

#ifdef DEBUG

   std::cout << "      value: " << sqrt(sqlength0) << " / " << sqrt(sqlength) << " = " << merit << std::endl;
#endif

   return pow(merit, m_power);
}

//===========================================================================

}

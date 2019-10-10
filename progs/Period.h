#ifndef LATMRG_PERIOD_PROGS_H
#define LATMRG_PERIOD_PROGS_H

extern std::ostream* out;

namespace LatMRG {
  template<typename Integ, typename Dbl> struct MaxPeriod {
    typedef Integ Int;
    GenType type = MRG;
    Int m_m, m_b, m_r;
    int m_k, m_e;

    NTL::vector<Int> m_a;
    Int m_mult;
    NTL::matrix<Int> m_A;

    DecompType decompm1;
    DecompType decompr;
    std::string filem1;
    std::string filer;

    int CheckPeriod()
    {
      MRGComponent<Int> mrg (m_m, m_k, decompm1, filem1.c_str(), decompr,
          filer.c_str());

      *out << "Period: A MRG Full Period Checker\n";
      *out << "==========================================================="
        "=====================\n\n";
      *out << "The generator with" << "\n";
      *out << "    m = " << m_m << " = " << m_b << "^" << m_e << "+"
        << m_r << "\n";
      *out << "    k = " << m_k << "\n";

      if (type == MRG) {
        *out << "    a = " << m_a << "\n";
        if (mrg.maxPeriod (m_a))
          *out << "has maximal period.\n";
        else
          *out << "does not have maximal period.\n";
      } else if (type == LCG) {
        *out << "    a = " << m_mult << "\n";
        if (mrg.maxPeriod (m_mult))
          *out << "can have maximal period 'm' if 'c' is choosen"
            " relatively prime to 'm'.\n";
        else
          *out << "cannot have maximal period.\n";
      } else if (type == MMRG) {
        *out << "    A = " << m_A << "\n";
        if (mrg.maxPeriod (m_A))
          *out << "has maximal period.\n";
        else
          *out << "does not have maximal period.\n";
      } else if (type == MWC) {
      }

      return 0;
    }
  };
}

#endif

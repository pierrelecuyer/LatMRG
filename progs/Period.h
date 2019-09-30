#ifndef LATMRG_PERIOD_PROGS_H
#define LATMRG_PERIOD_PROGS_H

namespace LatMRG {
  template<typename Integ, typename Dbl> struct MaxPeriod {
    typedef Integ Int;
    GenType type = MRG;
    Int m_m, m_b, m_r;
    int m_k, m_e;

    NTL::vector<Int> m_a;
    Int m_mult;
    NTL::vector<Int> m_A;

    DecompType decompm1;
    DecompType decompr;
    std::string filem1;
    std::string filer;

    int CheckPeriod()
    {
      MRGComponent<Int> mrg (m_m, m_k, decompm1, filem1.c_str(), decompr,
          filer.c_str());

      std::cout << "Period: A MRG Full Period Checker\n";
      std::cout << "==========================================================="
        "=====================\n\n";
      std::cout << "The generator with" << "\n";
      std::cout << "    m = " << m_m << " = " << m_b << "^" << m_e << "+"
        << m_r << "\n";
      std::cout << "    k = " << m_k << "\n";
      std::cout << "    a = " << m_a << "\n";

      if (mrg.maxPeriod (m_a))
        std::cout << "has maximal period.\n";
      else
        std::cout << "does not have maximal period.\n";

      return 0;
    }
  };
}

#endif

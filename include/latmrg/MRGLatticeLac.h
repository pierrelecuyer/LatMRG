#ifndef MRGLATTICELAC_H
#define MRGLATTICELAC_H

#include "latticetester/Const.h"
#include "latticetester/Lacunary.h"

#include "latmrg/Const.h"
#include "latmrg/MRGLattice.h"


namespace LatMRG {

  /**
   * This class implements lattice bases built from multiple recursive linear
   * congruential generators (see class <tt>MRGLattice</tt>) using *lacunary
   * indices*.
   *
   */
  template<typename Int, typename Dbl>
    class MRGLatticeLac:   public MRGLattice<Int, Dbl> {
      public:

        /**
         * Constructor with modulus of congruence \f$m\f$, order of the recurrence
         * \f$k\f$, multipliers \f$A\f$, maximal dimension `maxDim`, and lattice type
         * `latt`. Vector and matrix indices vary from 1 to `maxDim`. The length of
         * the basis vectors is computed with `norm`. The bases are built using the
         * *lacunary indices* `lac`.
         */
        MRGLatticeLac (const Int & m, const MVect & A, int maxDim, int k,
            BVect & lac, LatticeType latt,
            LatticeTester::NormType norm = LatticeTester::L2NORM);

        /**
         * Copy constructor. The maximal dimension of the new basis is set to
         * <tt>Lat</tt>’s current dimension.
         */
        MRGLatticeLac (const MRGLatticeLac<Int, Dbl> & Lat);

        /**
         * Assigns `Lat` to this object. The maximal dimension of this basis is
         * set to <tt>Lat</tt>’s current dimension.
         */
        MRGLatticeLac<Int, Dbl> & operator= (const MRGLatticeLac<Int, Dbl> & Lat);

        /**
         * Destructor.
         */
        virtual ~MRGLatticeLac();

        /**
         * Builds the basis of the MRG recurrence in dimension \f$d\f$ using
         * the lacunary indices.
         */
        void buildBasis (int d);

        /**
         * Increases the dimension of the basis by 1.
         */
        void incDim() { incDimBasis (this->getDim()); }

        /**
         * Returns the \f$j\f$-th lacunary index.
         */
        Int & getLac (int j);

        /**
         * Sets the lacunary indices for this lattice to `lat`.
         */
        void setLac (const LatticeTester::Lacunary<Int> & lat);

      protected:

        void initStates();

        /**
         * Increases the dimension of the basis by 1.
         */
        void incDimBasis (int);

        /**
         * The lacunary indices.
         */
        LatticeTester::Lacunary<Int> m_lac;
    }; // End class declaration

  //===========================================================================

  template<typename Int, typename Dbl>
    MRGLatticeLac<Int, Dbl>::MRGLatticeLac (const MRGLatticeLac<Int, Dbl> & lat): MRGLattice<Int, Dbl>::
                                                                        MRGLattice (lat.m_modulo, lat.m_aCoef, lat.getDim (), lat.getOrder (),
                                                                            lat.m_latType, lat.getNorm ()), m_lac (lat.m_lac)
  {
    this->m_lacunaryFlag = true;
    LatticeTester::MyExit (1, " MRGLatticeLac::  copy constructeur n'est pas terminé   ");
  }


  //=========================================================================

  template<typename Int, typename Dbl>
    MRGLatticeLac<Int, Dbl> & MRGLatticeLac<Int, Dbl>::operator= (const MRGLatticeLac<Int, Dbl> & lat)
    {
      if (this == &lat)
        return * this;
      LatticeTester::MyExit (1, " MRGLatticeLac::operator= n'est pas terminé   ");
      copyBasis (lat);
      return *this;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    MRGLatticeLac<Int, Dbl>::MRGLatticeLac (const Int & m, const MVect & a, int maxDim,
        int k, BVect & Lac, LatticeType lat,
        LatticeTester::NormType norm):
      MRGLattice<Int, Dbl>::MRGLattice (m, a, maxDim, k, lat, norm),
      m_lac (Lac, maxDim)
  {
    this->m_lacunaryFlag = true;
    this->m_sta.SetDims(1, 1);
    for (int i = 0; i < this->m_order; i++)
      this->m_aCoef[i] = a[i];
  }


  //============================================================================

  template<typename Int, typename Dbl>
    MRGLatticeLac<Int, Dbl>::~MRGLatticeLac ()
    {
    }


  //=========================================================================

  template<typename Int, typename Dbl>
    Int & MRGLatticeLac<Int, Dbl>::getLac (int j)
    {
      if (j <= this->m_lac.getSize () && j > 0)
        return this->m_lac.getLac (j);
      throw std::out_of_range ("MRGLatticeLac::getLac");
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void MRGLatticeLac<Int, Dbl>::setLac (const LatticeTester::Lacunary<Int> & lac)
    {
      this->m_lac = lac;
      this->m_lacunaryFlag = true;
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void MRGLatticeLac<Int, Dbl>::buildBasis (int d)
    {
      int ord = this->m_order;

      initStates ();
      int IMax = this->m_lac.getSize ();

      MVect b;
      b.SetLength (this->m_order+1);
      LatticeTester::Invert (this->m_aCoef, b, this->m_order);

      // b is the characteristic polynomial
      PolyPE<Int>::setM (this->m_modulo);
      PolyPE<Int>::setF (b);
      PolyPE<Int> pol;

      // Construction d'un systeme generateur modulo m.
      for (int k = 0; k < IMax; k++) {
        // pour chaque indice lacunaire
        NTL::conv (this->m_e, this->m_lac[k]);

        // x^m_e Mod f(x) Mod m
        pol.powerMod (this->m_e);
        pol.toVector (this->m_xi);

        for (int i = 0; i < this->m_order; i++) {
          this->m_wSI[i][k] = this->m_xi[i];
        }
      }


      /*
         cout << "m_order = " << m_order << endl;
         cout << "IMax = " << IMax << endl;
         cout << "m_lac = " << m_lac.toString() << endl;
         cout << "m_wSI = \n" << m_wSI << endl;
         cout << "m_vSI = \n" << m_vSI << endl;
         */


      //On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
      //les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
      //V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
      //avec m sur la diagonale).

      //Il serait possible de réserver m_wSI, m_vSI avec seulement IMax
      //lignes et colonnes quand order >> IMax. Mais il faudrait un nouveau
      //Triangularization qui pourrait nécessiter un appel pour chaque ligne,
      //mais qui sauverait beaucoup de mémoire.
      //Il n'est pas certain que cela en vaille la peine.
      LatticeTester::Triangularization <BMat> (this->m_wSI, this->m_vSI, ord, IMax, this->m_modulo);
      LatticeTester::CalcDual <BMat> (this->m_vSI, this->m_wSI, IMax, this->m_modulo);

      // Construire la base de dimension 1
      this->m_basis[0][0] = this->m_vSI[0][0];
      this->m_dualbasis[0][0] = this->m_wSI[0][0];
      this->setDim (1);

      this->setNegativeNorm ();
      this->setDualNegativeNorm ();

      for (int i = 1; i < d; i++)
        incDimBasis (IMax);

      // for debugging
      // trace("ESPION_2", i);
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void MRGLatticeLac<Int, Dbl>::incDimBasis (int IMax)
    {
      MRGLattice<Int, Dbl>::incDimLaBasis(IMax);
    }


  //===========================================================================

  template<typename Int, typename Dbl>
    void MRGLatticeLac<Int, Dbl>::initStates ()
    {
      /*
       * Initialise la matrice carrée Sta. La matrice Sta est d'ordre égal à l'ordre du
       * générateur. Elle contient un système de générateurs pour le groupe d'états
       * considérés.
       */
      NTL::clear (this->m_t2);

      if (this->m_latType == RECURRENT) {
        // check if a_k is relatively prime with m ==> m_t1 = 1
        this->m_t1 = NTL::GCD (this->m_aCoef[this->m_order], this->m_modulo);
        this->m_t1 = abs (this->m_t1);
        NTL::set (this->m_t2);
      }

      if (this->m_latType == FULL || this->m_latType == PRIMEPOWER || (this->m_t1 == this->m_t2)) {
        // m_sta is set to identity matrix
        this->calcLgVolDual2 (this->m_lgm2);

      } else {
        if (this->m_latType == ORBIT) {
          LatticeTester::MyExit (1, "case ORBIT ne fonctionne pas");
        } else if (this->m_latType == RECURRENT) {
          LatticeTester::MyExit (1, "case RECURRENT ne fonctionne pas");
        }
      }
    }

} // End namespace LatMRG
#endif

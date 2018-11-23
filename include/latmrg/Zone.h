#ifndef ZONE_H
#define ZONE_H

#include <iostream>
#include <cassert>

#include "latticetester/Util.h"

#include "latmrg/SeekConfig.h"

namespace LatMRG {

  /**
   * This class implements search zones in parameter space for the coefficients
   * of the recurrences defining generators or lattices.
   */
  template<typename Int>
    class Zone {
      public:

        /**
         * Possible zone number. <tt>ZONE1</tt> is the case where \f$b
         * < -\sqrt{m}\f$, <tt>ZONE3</tt> is the case where \f$b > \sqrt{m}\f$, and
         * <tt>ZONE2</tt> is the case where \f$ -\sqrt{m} \le b \le\sqrt{m}\f$.
         */
        enum ZoneType { ZONE1, ZONE2, ZONE3, NZONES };

        /**
         * Constructor.
         */
        Zone ();

        /**
         * Destructor.
         */
        ~Zone ();

        /**
         * Initializes the research zones for the multiplier \f$a_i\f$ of
         * <tt>comp</tt> which is the \f$s\f$-th component of the combined
         * generator. In the case of a random search, creates also the region
         * for this multiplier.
         */
        void init (const Component<Int> & comp, int s, int i);

        /**
         * Computes the lower bound <tt>inf</tt> of the intersection of the
         * zone with the search interval \f$[b,c]\f$. If <tt>z</tt> is equal to
         * <tt>ZONE1</tt> or <tt>ZONE3</tt>, the computed bound is such that
         * \f$\lfloor m/q\rfloor= c\f$ (approximately). If <tt>z = ZONE2</tt>,
         * the computed bound is \f$b\f$. In any case, the lower bound of the
         * zone is always \f$\le b\f$.
         */
        void calcInfBound (const Int & b, const Int & c, const Component<Int> & comp,
            int z, Int & inf);

        /**
         * Computes the upper bound <tt>sup</tt> of the intersection of the
         * zone with the search interval \f$[b,c]\f$. If <tt>z</tt> is equal to
         * <tt>ZONE1</tt> or <tt>ZONE3</tt>, the computed bound is such that
         * \f$\lfloor m/q\rfloor= b\f$ (approximately). If <tt>z = ZONE2</tt>,
         * the computed bound is \f$\min\{c, \sqrt{m}\}\f$.
         */
        void calcSupBound (const Int & b, const Int & c, const Component<Int> & comp,
            int z, Int & sup);

        /**
         * Returns the lower boundary of this region.
         */
        Int & getInf() { return inf; }

        /**
         * Returns the upper boundary of this region.
         */
        Int & getSup() { return sup; }
        Int & getSupMsH() { return supMsH; }

        void setSup(Int sup) { this->sup = sup;}
        void setInf(Int inf) { this->inf = inf;}

        /**
         * Returns the value of `frac`.
         */
        double getFrac() { return frac; }

        /**
         * Is <tt>true</tt> if and only if <tt>sup - inf \f$\le\f$ H</tt> (or
         * <tt>Hk</tt>).
         */
        bool smallF;

        /**
         * Returns the zone number for this region.
         */
        ZoneType getNo () { return No; }

        /**
         * Selects a random region and initializes its boundaries. The program
         * will search this region exhaustively.
         */
        void chooseBoundaries (const Component<Int> & comp, Zone<Int> *zone, long h);

        /**
         * Sets the values of the upper boundaries in the three zones based on
         * \f$m_j\f$ for each of the \f$J\f$ components generators. Also sets
         * the seed for the random number generator used in the random search.
         */
        static void initFrontieres (const SeekConfig<Int> & config);

        /**
         * Upper boundary of each zone.
         */
        static NTL::matrix<Int> Frontiere;

        /**
         * Is <tt>true</tt> if at upper boundary of each zone.
         */
        static const bool DivQ[NZONES];

        /**
         * List of zones.
         */
        Zone<Int> *nextZone;

        /**
         * Returns this zone as a string.
         */
        std::string toString();

      protected:

        /**
         * The zone number where this region is found.
         */
        ZoneType No;

        /**
         * Lower boundary defining the region.
         */
        Int inf;

        /**
         * Upper boundary defining the region.
         */
        Int sup;

        /**
         * The fraction of the acceptable values of the multiplier \f$a\f$
         * lying in this zone.
         */
        double frac;

        /**
         * When `small` is `true`, `supMsH` is equal to `sup+1-H` (or `sup+1-Hk`).
         */
        Int supMsH;

        /**
         * Work variables.
         */
        Int T1, T2;
    };

  //============================================================================

  template<typename Int>
  NTL::matrix<Int> Zone<Int>::Frontiere;

  // Si DivQ[i] = true, on cherche les a=m/q avec q entre les bornes.
  template<typename Int>
  const bool Zone<Int>::DivQ[] = {true, false, true}; // For ZONE1, ZONE2, ZONE3



  template<typename Int>
    Zone<Int>::Zone ()
    {
      // inf and sup at 1 is more logical since we will never really want zones at
      //0
      inf = 1;
      sup = 1;
      No = ZONE3;
      frac = 0.0;
      smallF = false;
      supMsH = 0;
      T1 = 0;
      T2 = 0;
      nextZone = NULL;
    }

  //===========================================================================

  template<typename Int>
    Zone<Int>::~Zone ()
    {
      if (nextZone != 0) {
        delete nextZone;
        nextZone = 0;
      }
      /*
         Zone<Int> *Z = nextZone;
         while (Z != 0) {
         Zone<Int> *p;
         p = Z->nextZone;
         delete Z;
         Z = p;
         }
         */
    }


  //===========================================================================

  template<typename Int>
    void Zone<Int>::initFrontieres (const SeekConfig<Int> & config)
    {
      LatticeTester::CreateMatr (Frontiere, config.J, NZONES);
      // Calcul des frontieres superieures des 3 zones
      for (int s = 0; s < config.J; s++) {
        Frontiere[s][ZONE1] = -(1 + config.compon[s].modulus.mRac);
        Frontiere[s][ZONE2] = config.compon[s].modulus.mRac;
        Frontiere[s][ZONE3] = config.compon[s].modulus.m - 1;
      }
      LatticeTester::SetSeed (config.seed);
    }

  //===========================================================================

  template<typename Int>
    void Zone<Int>::init (const Component<Int> & comp, int j, int i)
    {
      Int h(0), Eb, Ec;
      // Calcul des bornes pour a[j][i] ou q (ZONE 1 ou 3) dans chaque zone
      if (comp.searchMethod != EXHAUST) {
        if (i+1 == comp.k)
          h = comp.Hk;
        else
          h = comp.H;         // Taille des regions pour a[j][i]
      }
      Eb = comp.b[i];
      Ec = comp.c[i];

      if (comp.implemCond != APP_FACT) {
        // Pas de contrainte d'implantation; une seule grande zone
        inf = Eb;
        sup = Ec;                // Les bornes de la zone

        if (comp.searchMethod != EXHAUST) {
          // Ajuster taille des régions qui seront tirées au hasard
          frac = 1.0;
          T1 = sup - inf;

          smallF = T1 < h;
          if (!smallF)
            supMsH = sup - h + 1;
        }
        nextZone = 0;

      } else {  // La contrainte AppFact est appliquee; Max de 3 zones
        Zone<Int> *Z = 0, *pZ = 0;
        Int Tot;
        Tot = 0;
        bool isFirst = true;
        int k = ZONE1;
        // Trouver la premiere zone k intersectant [Eb, Ec]
        while (k < NZONES - 1 && Frontiere[j][k] < Eb)
          ++k;

        // Initialisation de la zone k
        while (Eb <= Ec) {
          if (isFirst) {
            isFirst = false;
            Z = this;
          } else {
            // pZ precedes Z in the list of zones
            pZ->nextZone = Z = new Zone<Int> ();
            // setFrontieres (comp.modulus);
            if (comp.searchMethod != EXHAUST)
              ; // comp.ApproxTotGen = true;
          }
          pZ = Z;
          Z->No = (ZoneType) k;
          Z->calcInfBound (Eb, Ec, comp, k, Z->inf); // Bornes de la zone
          Z->calcSupBound (Eb, Ec, comp, k, Z->sup);
          T1 = Z->sup - Z->inf;        // Taille de la zone, -1
          Tot += T1;
          ++Tot;
          if (comp.searchMethod != EXHAUST) {
            // Ajuster taille des régions qui seront tirées au hasard
            Z->smallF = T1 < h;
            if (Z->smallF) {
              std::cout << "** WARNING ** Multiplier " << i <<
                " zone " << k << " : sup - inf + 1 < H(k). " << std::endl;
            } else {
              Z->supMsH = Z->sup - h + 1;
            }
          }
          Eb = Frontiere[j][k] + 1; // passer a la zone suivante
          ++k;
        }                        // while
        Z->nextZone = 0;

        // Calculer la proportion frac que represente chaque zone
        Z = this;
        double LRE, TotLR;
        while (Z != 0) {
          T1 = Z->sup - Z->inf + 1;
          NTL::conv (LRE, T1);
          NTL::conv (TotLR, Tot);
          Z->frac = LRE / TotLR;
          Z = Z->nextZone;
        }
      }
    }


  //===========================================================================

  template<typename Int>
    void Zone<Int>::calcInfBound (const Int & b, const Int & c,
        const Component<Int> & comp, int No, Int & inf)
    {
      /* Calcule la borne inferieure de l'intersection d'une zone avec [b..c].
         Si No = 1 ou 3, la borne retournee sera (approx) q t.q. [m/q] = c
         Si No = 2, retourne b. A l'appel, b est toujours >= borne inf. de
         la zone. */

      switch (No) {
        case ZONE1:
          if (c < comp.modulus.mRacNeg) {
            T1 = comp.modulus.m - 1;
            T2 = c + 1;
            LatticeTester::Divide (T1, T2, inf, T1);
            if (!NTL::IsZero (T1))
              inf += 1;
            if (inf < (comp.modulus.mRacNeg))
              inf = comp.modulus.mRacNeg;
          } else
            inf = comp.modulus.mRacNeg;
          break;

        case ZONE2:                      // -sqrt(m)-1 < b <= sqrt(m);
          inf = b;
          break;

        case ZONE3:
          /* b > sqrt(m) On retourne le a minimal tel que b <= floor(m/a) <= c.
Note: floor(m/Borne) <= c, m/Borne < c+1, Borne > m/(c+1). On fait:
Borne = (m DIV (c+1)) + 1. */
          T1 = c + 1;
          LatticeTester::Quotient (comp.modulus.m, T1, T2);
          inf = T2 + 1;
          if (inf == 1)               // si inf = 1, alors a = m
            inf = 2;
          break;

        default:
          std::cout << "Zone::calcInfBound:   impossible case" << std::endl;
          assert (0);
      }
    }


  //===========================================================================

  template<typename Int>
    void Zone<Int>::calcSupBound (const Int & b, const Int & c,
        const Component<Int> & comp, int No, Int & sup)
    {
      switch (No) {
        case ZONE1:                      // b <= -sqrt(m)-1
          LatticeTester::Quotient (comp.modulus.m, b, sup);
          if (sup == -1)               // si sup = -1, alors a = -m
            sup = -2;
          break;

        case ZONE2:                      // -sqrt(m)-1 < b <= sqrt(m)
          if (comp.modulus.mRac < c)
            sup = comp.modulus.mRac;
          else
            sup = c;
          break;

        case ZONE3:                      /* b > sqrt(m) On retourne le a maximal tel
                                            que b <= floor(m/a) <= c. */
          LatticeTester::Quotient (comp.modulus.m, b, sup);
          if (comp.modulus.mRac < sup)
            sup = comp.modulus.mRac;
          break;

        default:
          std::cout << "Zone::calcSupBound:   impossible case" << std::endl;
          assert (0);
      }
    }


  //===========================================================================

  template<typename Int>
    std::string Zone<Int>::toString()
    {
      std::ostringstream sortie;
      sortie << "No = " << No << std::endl
        << "inf = " << inf << std::endl
        << "sup = " << sup << std::endl
        << "frac = " << frac << std::endl
        << "supMsH = " << supMsH << std::endl
        << "smallF = " << smallF << std::endl
        << "DivQ = { " << DivQ[0] << ", " << DivQ[1] << ", " <<
        DivQ[2] << " }" << std::endl
        << "Frontiere = { " << Frontiere[0][0] << ", " << Frontiere[0][1] << ", " <<
        Frontiere[0][2] << " }" << std::endl
        << nextZone << std::endl << std::endl;

      return sortie.str ();
    }


  //===========================================================================
  /*
     void Zone<Int>::toRegion (Zone<Int> *Z)
     {
     inf = Z->getInf();
     sup = Z->getSup();
     No = Z->getNo ();
     frac = Z->getFrac();
     smallF = Z->smallF;
     supMsH = Z->getSupMsH();
     cond = Z->cond;
     nextZone = Z->nextZone;
     Frontiere[0] = Z->Frontiere[0];
     Frontiere[1] = Z->Frontiere[1];
     Frontiere[2] = Z->Frontiere[2];
     }
     */

  //===========================================================================

  template<typename Int>
    void Zone<Int>::chooseBoundaries (const Component<Int> & comp, Zone<Int> *zone, long h)
    {
      /* Choisir une region au hasard et en initialiser les bornes.
       * On fouillera ensuite completement cette region.
       * Utilise avec la methode de recherche RANDOM seulement. */

      Zone<Int> *Z = zone;
      double p, u;
      Zone<Int> *R = this;

      //for (int i = 0; i < comp.k; ++i) {
      if (Z->nextZone != 0) {
        // Choisir d'abord l'une des zones selon les probabilites Z->frac
        u = LatticeTester::RandU01();
        p = Z->getFrac();
        while ((u > p) && (Z != 0)) {
          Z = Z->nextZone;
          p += Z->getFrac();
        }
      }


      // Z pointe a la zone choisie
      if (Z->smallF) {
        // On prends toute la zone
        R->setInf(Z->getInf());
        R->setSup(Z->getSup());
      } else {
        // On prends un intervalle de taille h, au hasard, dans la zone
        /*
           if (i == comp.k)
           h = comp.Hk - 1;
           else
           h = comp.H - 1;*/
        T1 = Z->getSupMsH() - Z->getInf() + 1;
        NTL::conv (p, T1);
        p *= LatticeTester::RandU01();
        NTL::conv (T1, p);
        T1 = Z->getInf() + T1;
        R->setInf(T1);
        R->setSup(T1 + h);
      }
      R->No = Z->getNo();
      //}
    }
}
#endif

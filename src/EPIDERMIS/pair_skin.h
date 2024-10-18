/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(skin,PairSkin);
// clang-format on
#else

#ifndef LMP_PAIR_SKIN_H
#define LMP_PAIR_SKIN_H

#include "pair.h"

namespace LAMMPS_NS {

    class PairSkin : public Pair {
     public:
      PairSkin(class LAMMPS *);
      ~PairSkin() override;

      void compute(int, int) override;
      void settings(int, char **) override;
      void coeff(int, char **) override;
      double init_one(int, int) override;

      void write_restart(FILE *) override;
      void read_restart(FILE *) override;
      void write_data(FILE *) override;
      void write_data_all(FILE *) override;

     protected:
      double cut_global;
      double **k_n;    // repulsion strength

      virtual void allocate();
    };

}    // namespace LAMMPS_NS

#endif
#endif

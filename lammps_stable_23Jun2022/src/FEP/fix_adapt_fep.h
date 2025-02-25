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

#ifdef FIX_CLASS
// clang-format off
FixStyle(adapt/fep,FixAdaptFEP);
// clang-format on
#else

#ifndef LMP_FIX_ADAPT_FEP_H
#define LMP_FIX_ADAPT_FEP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdaptFEP : public Fix {
 public:
  int diamflag;    // 1 if atom diameters will vary, for AtomVecGranular
  int chgflag;

  FixAdaptFEP(class LAMMPS *, int, char **);
  ~FixAdaptFEP() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void post_run() override;
  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  void set_arrays(int) override;

 private:
  int nadapt, resetflag, scaleflag, afterflag;
  int anypair;
  int nlevels_respa;
  char *id_fix_diam, *id_fix_chg;
  class FixStore *fix_diam, *fix_chg;

  struct Adapt {
    int which, ivar;
    char *var;
    char *pstyle, *pparam;
    int ilo, ihi, jlo, jhi;
    int pdim;
    double *scalar, scalar_orig;
    double **array, **array_orig;
    int aparam;
  };

  Adapt *adapt;
  double *kspace_scale;

  void change_settings();
  void restore_settings();
};

}    // namespace LAMMPS_NS

#endif
#endif

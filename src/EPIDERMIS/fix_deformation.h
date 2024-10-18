/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(skin/deformation,FixDeformation)

#else

#ifndef LMP_FIX_DEFORMATION_H
#define LMP_FIX_DEFORMATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeformation: public Fix {
 public:
  FixDeformation(class LAMMPS *, int, char **);
  virtual ~FixDeformation();

  int setmask();
  void biology_nufeb();
  void compute(int);

 protected:

  double alpha_f, beta_f, p_f;

  double *diffstate;        // per-type list of differentiated cell shape
  int *diffgroup;
  int group_flag;
  int *groupbits;
  int ngroups;

  class FixPropertyDiffstate *fix_diffstate;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

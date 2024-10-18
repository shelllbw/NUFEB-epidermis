/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04- 94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(skin/differentiation/ta,FixDifferentiationTA)

#else

#ifndef LMP_FIX_CELL_DIFFERENTIATION_TA_H
#define LMP_FIX_CELL_DIFFERENTIATION_TA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDifferentiationTA: public Fix {
 public:
  FixDifferentiationTA(class LAMMPS *, int, char **);
  virtual ~FixDifferentiationTA();

  virtual void init();
  void post_constructor();
  int setmask();
  void biology_nufeb();
  void compute();
  double get_zhi(int);

 protected:
  double ave_time, sd;     // average differentiation time with sd

  int icyto;
  int type_spin;    // atom type of spinous cell
  int mask_spin;    // Spinous group mask
  int seed;
  char *group_id;
  double cyto_affinity;

  class RanPark *random;
  class FixPropertyCycletime *fix_ct;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

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

FixStyle(skin/division/stem,FixDivideStem)

#else

#ifndef LMP_FIX_DIVIDE_STEM_H
#define LMP_FIX_DIVIDE_STEM_H

#include "fix_divide.h"

namespace LAMMPS_NS {

class FixDivideStem : public FixDivide {
 public:

  FixDivideStem(class LAMMPS *, int, char **);
  virtual ~FixDivideStem();
  virtual void compute();
  virtual void init();
  
 protected:
  int type_ta;    // atom type of TA cell
  int mask_ta;    // TA group mask
  int seed;
  char *group_id;

  double pa1, pa2, pb1, pb2;

  class RanPark *random;
  class FixPropertyCycletime *fix_ct;
};
}

#endif
#endif

/* ERROR/WARNING messages:
*/

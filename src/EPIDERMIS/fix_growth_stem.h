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

FixStyle(nufeb/growth/stem,FixGrowthStem)

#else

#ifndef LMP_FIX_GROWTH_STEM_H
#define LMP_FIX_GROWTH_STEM_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthStem: public FixGrowth {
 public:
  FixGrowthStem(class LAMMPS *, int, char **);
  virtual ~FixGrowthStem() {}

  virtual void update_atoms();
  virtual void update_cells() {}

 protected:
  int isub;
  double growth;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

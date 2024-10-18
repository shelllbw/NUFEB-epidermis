/* ----------------------------------------------------------------------
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

FixStyle(skin/adhesion,FixAdhesionSkin)

#else

#ifndef LMP_FIX_ADHESION_SKIN_H
#define LMP_FIX_ADHESION_SKIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdhesionSkin : public Fix {
 public:
  class NeighList *list;

  FixAdhesionSkin(class LAMMPS *, int, char **);
  virtual ~FixAdhesionSkin();
  void allocate();
  void init();
  void init_list(int, class NeighList *);
  int modify_param(int, char **);
  int setmask();
  virtual void post_force(int);

 protected:
  double **af; //adhesion factor
  double zeta;
  double r;
  int allocated;

  void compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

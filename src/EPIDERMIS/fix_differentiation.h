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

FixStyle(epidermis/differentiation,FixDifferentiation)

#else

#ifndef LMP_FIX_CELL_DIFFERENTIATION_H
#define LMP_FIX_CELL_DIFFERENTIATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDifferentiation: public Fix {
 public:
  FixDifferentiation(class LAMMPS *, int, char **);
  virtual ~FixDifferentiation();

  int setmask();
  void biology_nufeb();
  void compute();

 protected:

  double *diffshape;        // per-type list of differentiated cell shape
  int *diffgroup;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

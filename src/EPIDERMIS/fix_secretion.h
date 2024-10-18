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

FixStyle(skin/secretion,FixSecretion)

#else

#ifndef LMP_FIX_SECRETION_H
#define LMP_FIX_SECRETION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSecretion: public Fix {
 public:
  FixSecretion(class LAMMPS *, int, char **);
  virtual ~FixSecretion() {}

  int setmask();
  void biology_nufeb();
  void secrete_gf();
  void secrete_cyto();

 protected:
  int isub, icyto;
  int tstem, tta, flag;
  double yield, decay;
  double cyto_max;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

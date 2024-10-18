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

#ifdef COMPUTE_CLASS

ComputeStyle(skin/base_height,ComputeBaseHeight)

#else

#ifndef LMP_COMPUTE_BASE_HEIGHT_H
#define LMP_COMPUTE_BASE_HEIGHT_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeBaseHeight : public Compute {
 public:
  ComputeBaseHeight(class LAMMPS *, int, char **);
  virtual ~ComputeBaseHeight() {}
  virtual void init() {}
  virtual double compute_scalar();

private:
  bigint nnormal;
  int istem, ita;
  double h_max, h_min;
  double tau;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

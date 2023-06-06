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

FixStyle(nve/bacillus/limit,FixNVEBacillusLimit)

#else

#ifndef LMP_FIX_BACILLUS_NVE_LIMIT_H
#define LMP_FIX_BACILLUS_NVE_LIMIT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVEBacillusLimit : public Fix {
 public:
  FixNVEBacillusLimit(class LAMMPS *, int, char **);
  virtual ~FixNVEBacillusLimit() {}
  int setmask();
  void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  void reset_dt();
  virtual double compute_scalar();

 protected:
  double dtv,dtf;
  double dtq;
  double *step_respa;
  int ncount;
  double xlimit,vlimitsq;

  class AtomVecBacillus *avec;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Should not use fix nve/limit with fix shake or fix rattle

This will lead to invalid constraint forces in the SHAKE/RATTLE
computation.

*/
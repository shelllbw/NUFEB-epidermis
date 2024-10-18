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

FixStyle(skin/update_base,FixUpdateBase)

#else

#ifndef LMP_FIX_UPDATE_BASE_H
#define LMP_FIX_UPDATE_BASE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixUpdateBase: public Fix {
 public:
  FixUpdateBase(class LAMMPS *, int, char **);
  virtual ~FixUpdateBase();

  virtual void biology_nufeb();
  int setmask();

  void add_lattice();
  void loop_lattice(int);
  void delete_base();

 private:
  int ntype, nbasis;
  int *basistype;
  int mask_base;

  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  int ilo, ihi, jlo, jhi, klo, khi;

  int *dlist;

  int nlatt;             // number of owned lattice sites
  int nlatt_overflow;    // 1 if local nlatt exceeds a 32-bit int

  class Region *region;

  int triclinic;
  double sublo[3], subhi[3];    // epsilon-extended proc sub-box for adding atoms

  int vartest(double *);    // evaluate a variable with new atom position
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

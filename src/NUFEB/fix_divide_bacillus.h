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

FixStyle(nufeb/division/bacillus,FixDivideBacillus)

#else

#ifndef LMP_FIX_DIVIDE_BACILLUS_H
#define LMP_FIX_DIVIDE_BACILLUS_H

#include "fix_divide.h"

namespace LAMMPS_NS {

class FixDivideBacillus : public FixDivide {
 public:
  
  FixDivideBacillus(class LAMMPS *, int, char **);
  virtual ~FixDivideBacillus();
  virtual void compute();
  
 private:
  int seed;
  int conserveflag;
  double var;

  class RanPark *random;
  class AtomVecBacillus *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/

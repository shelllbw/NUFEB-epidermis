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

FixStyle(epidermis/division/stem,FixDivideStem)

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
  
 protected:
    int type_ta;    // atom type of TA cell
    int ita;       // list of groups contains TA cells
    int seed;

  void get_location(int, double*, double*, double, double);

  class RanPark *random;
  class FixPropertyCycletime *fix_ct;
};
}

#endif
#endif

/* ERROR/WARNING messages:
*/

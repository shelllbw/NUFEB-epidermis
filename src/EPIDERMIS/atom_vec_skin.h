/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(skin,AtomVecSkin);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_SKIN_H
#define LMP_ATOM_VEC_SKIN_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecSkin : public AtomVec {
 public:
  AtomVecSkin(class LAMMPS *);
  void process_args(int, char **);
  void init();

  void grow_pointers();
  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

 private:
  double *rmass, *biomass;    // wet and dry mass
  double **shape;             // 3 radii in xyz
  double *alpha;              // deformation rate

  int skin_flag;
  double rmass_one;
  double shape_a, shape_b, shape_c;
  double alpha_one;

  int shapevary;
};

}    // namespace LAMMPS_NS

#endif
#endif

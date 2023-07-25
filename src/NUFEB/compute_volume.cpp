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

#include "compute_volume.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define FOURTHIRDSPI 4.18879020478

/* ---------------------------------------------------------------------- */

ComputeVolume::ComputeVolume(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute nufeb/volume command");
  avec = nullptr;
  avec = (AtomVecBacillus *) atom->style_match("bacillus");

  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeVolume::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  int *mask = atom->mask;

  scalar = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    if ((mask[i] & groupbit)) {
      if (atom->bacillus_flag) {
        double r = atom->radius[i];
        int ibonus = atom->bacillus[i];
        AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

        double vsphere = FOURTHIRDSPI * r * r * r;
        double vcylinder = MY_PI * r * r * bonus->length;
        scalar += vsphere + vcylinder;
      } else if (atom->coccus_flag) {
        double r = atom->radius[i];
        scalar += FOURTHIRDSPI * r * r * r;
      } else if (atom->skin_flag) {
        scalar += FOURTHIRDSPI * atom->shape[i][0] * atom->shape[i][1] * atom->shape[i][2];
      } else {
        error->all(FLERR,"Unsupport atom types");
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  return scalar;
}

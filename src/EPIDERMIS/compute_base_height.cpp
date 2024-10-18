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

#include "compute_base_height.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeBaseHeight::ComputeBaseHeight(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute skin/base_height command");

  tau = utils::numeric(FLERR,arg[3],true,lmp);
  if (tau < 0)
    error->all(FLERR, "Illegal compute skin/base_height command: tau");

  istem = group->find(arg[4]);
  ita = group->find(arg[5]);

  h_max = 1.26e-4;
  h_min = 0.4e-4;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "h_min") == 0) {
      h_max = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "h_max") == 0) {
      h_min = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/division/bacillus command");
    }
  }

  scalar_flag = 1;
  scalar = h_min;

  nnormal = group->count(istem) + group->count(ita);
}

/* ---------------------------------------------------------------------- */

double ComputeBaseHeight::compute_scalar()
{
  double r, h;

  invoked_scalar = update->ntimestep;

  r = (double)nnormal / ((double)group->count(istem) + (double)group->count(ita));

  if (r > 1) r = 1;
  else if (r < 0) r = 0;

  h = h_max - (h_max - h_min) * r;
  scalar = scalar + (h - scalar) / tau * update->dt;

  return scalar;
}

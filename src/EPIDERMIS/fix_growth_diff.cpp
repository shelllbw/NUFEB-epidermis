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

#include <cstring>
#include "atom.h"
#include "error.h"

#include "fix_growth_diff.h"
#include "grid.h"
#include "group.h"
#include "modify.h"
#include "grid_masks.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGrowthDiff::FixGrowthDiff(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix epidermis/growth/diff command");

  ical = -1;
  deform = 0.0;
  yield = 1.0;

  ical = grid->find(arg[3]);
  if (ical < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "deform") == 0) {
      deform = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0 ) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix epidermis/growth/diff command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthDiff::update_atoms()
{
  double **conc = grid->conc;
  // use growth variable to store specific deformation rate
  double ***growth = grid->growth;
  int nlocal = atom->nlocal;
  double **shape = atom->shape;
  double *alpha = atom->alpha;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double tmp;
      int cell = grid->cell(atom->x[i]);
      // update cell shape based on flattening rate
      tmp = pow(conc[ical][cell] / (1e-6), 2);
      alpha[i] = (1 + deform * tmp) / (1 + tmp);

      printf("appha = %e %e \n", alpha[i], conc[ical][cell]);
      shape[i][0] *= alpha[i];
      shape[i][1] *= alpha[i];
      shape[i][2] /= alpha[i] * alpha[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthDiff::update_cells ()
{
  double **reac = grid->reac;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      // calcium secretion
      int cell = grid->cell(atom->x[i]);
      reac[ical][cell] += yield;
    }
  }
}
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

#include "fix_growth_aob.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowthAOB::FixGrowthAOB(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/growth/aob command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/aob requires grid_style nufeb/chemostat");

  inh4 = -1;
  io2 = -1;
  ino2 = -1;

  nh4_affinity = 0.0;
  o2_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;

  inh4 = grid->find(arg[3]);
  if (inh4 < 0)
    error->all(FLERR, "Can't find substrate name");
  nh4_affinity = utils::numeric(FLERR,arg[4],true,lmp);

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate name");
  o2_affinity = utils::numeric(FLERR,arg[6],true,lmp);

  ino2 = grid->find(arg[7]);
  if (ino2 < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maintain") == 0) {
      maintain = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/aob command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthAOB::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      double tmp1 = growth * conc[inh4][i] / (nh4_affinity + conc[inh4][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
      double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

      reac[inh4][i] -= 1 / yield * tmp1 * dens[igroup][i];
      reac[io2][i] -= (3.42 - yield) / yield * tmp1 * dens[igroup][i] + tmp2 * dens[igroup][i];
      reac[ino2][i] += 1 / yield * tmp1 * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthAOB::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = growth * conc[inh4][i] / (nh4_affinity + conc[inh4][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
    double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    grid->growth[igroup][i][0] = tmp1 - tmp2 - decay;
  }

  update_atoms_coccus();
}

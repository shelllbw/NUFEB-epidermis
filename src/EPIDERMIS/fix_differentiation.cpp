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

#include "fix_differentiation.h"
#include "atom.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDifferentiation::FixDifferentiation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/differentiation requires skin atom style");

  if (narg < 3)
    error->all(FLERR, "Illegal fix epidermis/differentiation command");

  diffshape = nullptr;
  diffgroup = nullptr;

  memory->grow(diffshape,atom->ntypes+1,"fix epidermis/differentiation:diffshape");
  memory->grow(diffgroup,atom->ntypes+1,"fix epidermis/differentiation:diffgroup");

  for (int i = 1; i <= atom->ntypes; i++) {
    diffshape[i] = MAXFLOAT;
  }

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "adiff") == 0) {
      int igroup = group->find(arg[iarg+1]);
      int type = utils::inumeric(FLERR,arg[iarg+2],true,lmp);
      double shapez = utils::numeric(FLERR,arg[iarg+3],true,lmp);

      if (igroup == -1)
        error->all(FLERR, "Can't find group in fix epidermis/differentiation command");

      diffshape[type] = shapez;
      diffgroup[type] = 1 | group->bitmask[igroup];
      iarg += 4;
    } else {
      error->all(FLERR, "Illegal fix epidermis/differentiation command");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixDifferentiation::~FixDifferentiation()
{
  memory->destroy(diffshape);
  memory->destroy(diffgroup);
}

/* ---------------------------------------------------------------------- */

int FixDifferentiation::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDifferentiation::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixDifferentiation::compute()
{
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int *mask = atom->mask;

  double **shape = atom->shape;
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      // reassign atom type based on shapez criterion
      for (int j = 1; j <= ntypes; j++) {
        if (shape[i][2] < diffshape[j]) {
          type[i] = j;
          mask[i] = diffgroup[j];
        }
      }
      printf("type[i] = %i \n", type[i]);
    }
  }
}
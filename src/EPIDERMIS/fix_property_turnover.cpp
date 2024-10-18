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

#include "fix_property_turnover.h"

#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "grid.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyTurnover::FixPropertyTurnover(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/turnover command");

  create_attribute = 1;
  // use vprop if size_peratom_cols = 1
  size_peratom_cols = 0;
  grow_arrays(atom->nmax);

  type_flag = 0;
  ntypes = 0;
  scalar_flag = 1;
  turnover_flag = 0;
  time = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "types") == 0) {
      type_flag = 1;
      ntypes = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      memory->create(atypes, ntypes, "fix skin/turnover:groupbits");

      for (int i = 0; i < ntypes; i++) {
        atypes[i] = utils::inumeric(FLERR, arg[iarg + 2 + i], true, lmp);
      }

      iarg = iarg + 2 + ntypes;
    } else if (strcmp(arg[iarg], "etype") == 0) {
      etype = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg = iarg + 2;
    } else {
      error->all(FLERR, "Illegal fix skin/property/turnover command");
    }
  }

  for (int i = 0; i < atom->nlocal; i++) {
    vprop[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */
FixPropertyTurnover::~FixPropertyTurnover()
{
  if (type_flag)
    memory->destroy(atypes);
}

/* ---------------------------------------------------------------------- */

int FixPropertyTurnover::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropertyTurnover::biology_nufeb()
{
  turnover_flag = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    for (int j = 0; j < ntypes; j++) {
      if (atom->type[i] == atypes[j]) {
        if (vprop[i] > 0.5) {
          turnover_flag = 1;
          break;
        }
      }
    }
  }

  if (!turnover_flag) {
    printf("turnover = %e type = %i \n", update->ntimestep * update->dt - time, atypes[0]);
    for (int i = 0; i < atom->nlocal; i++) {
      vprop[i] = 0.0;
    }
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = 0; j < ntypes; j++) {
        if (atom->type[i] == atypes[j]) {
          vprop[i] = 1.0;
        }
      }
    }
    time = update->ntimestep * update->dt;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyTurnover::set_arrays(int j)
{
  vprop[j] = 0.0;
}

/* ----------------------------------------------------------------------
   update array values of two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyTurnover::update_arrays(int i, int j)
{
  vprop[j] = 0.0;
  vprop[i] = 0.0;
}



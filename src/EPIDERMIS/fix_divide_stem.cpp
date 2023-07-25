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

#include "fix_divide_stem.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "compute.h"
#include "random_park.h"
#include "modify.h"
#include "domain.h"
#include "grid.h"
#include "group.h"
#include "fix_property_cycletime.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDivideStem::FixDivideStem(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/division/stem requires skin atom style");

  if (narg < 6)
    error->all(FLERR, "Illegal fix epidermis/division/stem command");

  type_ta = utils::inumeric(FLERR,arg[3],true,lmp);
  ita = group->find(arg[4]);
  if (ita < 0)
    error->all(FLERR, "Can't find TA group in fix epidermis/division/stem command");
  seed = utils::inumeric(FLERR,arg[5],true,lmp);

  fix_ct = nullptr;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"epidermis/property/cycletime") == 0) {
      fix_ct = (FixPropertyCycletime *) modify->fix[i];
    }
  if (fix_ct == nullptr)
    error->all(FLERR,"fix epidermis/division/stem requires fix epidermis/property/cycletime");

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

/* ---------------------------------------------------------------------- */

FixDivideStem::~FixDivideStem()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixDivideStem::compute()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell =grid->cell(x[i]);
      double growth = grid->growth[igroup][cell][0];
      double div_time = log2(growth);

      // trigger cell division if current cell cycle time
      // is greater than calculated doubling time
      if (fix_ct->aprop[i][0] > div_time) {

      }
    }
  }
}

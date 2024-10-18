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

#include "fix_differentiation_ta.h"
#include "atom.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "grid.h"
#include "modify.h"
#include "update.h"
#include "comm.h"
#include "random_park.h"
#include "fix_property_cycletime.h"
#include "fix_property_diffstate.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDifferentiationTA::FixDifferentiationTA(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/differentiation requires skin atom style");

  if (narg < 9)
    error->all(FLERR,  "Illegal fix epidermis/differentiation command");

  type_spin = utils::inumeric(FLERR,arg[3],true,lmp);
  mask_spin = group->find(arg[4]);
  mask_spin = 1 | group->bitmask[mask_spin];


  if (mask_spin < 0)
    error->all(FLERR, "Can't find TA group in fix epidermis/division/stem command");

  ave_time = utils::numeric(FLERR,arg[5],true,lmp);
  sd = utils::numeric(FLERR,arg[6],true,lmp);

  icyto = grid->find(arg[7]);
  if (icyto < 0)
    error->all(FLERR, "Can't find cytokine name");

  cyto_affinity = 0.0;
  cyto_affinity = utils::numeric (FLERR,arg[8],true,lmp);
  if (cyto_affinity <= 0)
    error->all(FLERR, "cytokine affinity must be greater than zero");

  seed = utils::inumeric(FLERR,arg[9],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  group_id = new char[strlen(arg[1])+1];
  strcpy(group_id, arg[1]);
}

/* ---------------------------------------------------------------------- */

FixDifferentiationTA::~FixDifferentiationTA()
{
  delete random;
  delete [] group_id;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */

void FixDifferentiationTA::post_constructor() {
  // create fix nufeb/property/cycletime
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"diff_ta_ct";
  fixarg[1] = group_id;
  fixarg[2] = (char *)"nufeb/property/cycletime";

  modify->add_fix(3, fixarg, 1);
  fix_ct = (FixPropertyCycletime *)modify->fix[modify->nfix-1];

  delete [] fixarg;
}

/* ---------------------------------------------------------------------- */

void FixDifferentiationTA::init()
{
//  for (int i = 0; i < atom->nlocal; i++) {
//    if (atom->mask[i] & groupbit) {
//      double prob = random->uniform();
//      double div_time = log(2) / 0.00001284;
//      fix_ct->aprop[i][0] = prob * div_time;
//      fix_ct->aprop[i][1] = prob * 216000;
//    }
//  }
}

/* ---------------------------------------------------------------------- */

int FixDifferentiationTA::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDifferentiationTA::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixDifferentiationTA::compute()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  double *bulk = grid->bulk;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      // double t = (ave_time + random->uniform()*sd - ave_time * 0.9 * (bulk[icyto] / (bulk[icyto] + 10)));
      double t = ave_time + random->gaussian() * sd;
      double t_cyto = t * 0.75 * bulk[icyto] / (bulk[icyto] + cyto_affinity);
      double lifetime = fix_ct->aprop[i][1];

      if (lifetime > t - t_cyto) {
        type[i] = type_spin;
        mask[i] = mask_spin;
        atom->x[i][2] = get_zhi(i);
      }
    }
  }
}

double FixDifferentiationTA::get_zhi(int i)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double zhi = 0.0;

  for (int j = 0; j < nlocal; j++) {
    if (atom->mask[j] & groupbit) {
      if (x[j][0] > x[i][0]-1e-5 && x[j][0] < x[i][0]+1e-5 &&
      x[j][1] > x[i][1]-1e-5 && x[j][1] < x[i][1]+1e-5)
        if (x[j][2] > zhi)
          zhi = x[j][2] + 5e-6;
    }
  }

  return zhi;
}

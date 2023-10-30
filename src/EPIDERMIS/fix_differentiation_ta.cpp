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
#include "modify.h"
#include "update.h"
#include "random_park.h"
#include "fix_property_cycletime.h"
#include "fix_property_generation.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDifferentiationTA::FixDifferentiationTA(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/differentiation requires skin atom style");

  if (narg < 6)
    error->all(FLERR,  "Illegal fix epidermis/differentiation command");

  type_spin = utils::inumeric(FLERR,arg[3],true,lmp);
  mask_spin = group->find(arg[4]);
  mask_spin = 1 | group->bitmask[mask_spin];

  if (mask_spin < 0)
    error->all(FLERR, "Can't find TA group in fix epidermis/division/stem command");

  ave_time = utils::numeric(FLERR,arg[5],true,lmp);
  sd = utils::numeric(FLERR,arg[6],true,lmp);

  seed = utils::inumeric(FLERR,arg[7],true,lmp);

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
   initialization before run
------------------------------------------------------------------------- */

void FixDifferentiationTA::init()
{
  // create fix nufeb/property/cycletime
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"diff_ta_ct";
  fixarg[1] = group_id;
  fixarg[2] = (char *)"nufeb/property/cycletime";

  modify->add_fix(3, fixarg, 1);
  delete [] fixarg;
  fix_ct = (FixPropertyCycletime *)modify->fix[modify->nfix-1];

  // create fix nufeb/property/generation
  char **fixarg2 = new char*[3];
  fixarg2[0] = (char *)"diff_ta_gen";
  fixarg2[1] = group_id;
  fixarg2[2] = (char *)"nufeb/property/generation";

  modify->add_fix(3, fixarg2, 1);
  delete [] fixarg2;
  fix_gen = (FixPropertyGeneration *)modify->fix[modify->nfix-1];
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

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double t = ave_time+ random->gaussian()*sd;
      double lifetime = fix_ct->aprop[i][0]+fix_ct->aprop[i][1]*fix_gen->vprop[i];
      if(atom->tag[i]==11500)printf("i=%i ct=%e  \n", i,lifetime);
      if (lifetime > t) {
        type[i] = type_spin;
        mask[i] = mask_spin;
      }
    }
  }
}

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

#include "fix_divide_ta.h"

#include "math_const.h"
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
using namespace MathConst;

#define DELTA 1.05

/* ---------------------------------------------------------------------- */

FixDivideTA::FixDivideTA(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix skin/division/ta requires skin atom style");

  if (narg < 4)
    error->all(FLERR, "Illegal fix skin/division/ta command");

  seed = utils::inumeric(FLERR,arg[3],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  group_id = new char[strlen(arg[1])+1];
  strcpy(group_id, arg[1]);
}

/* ---------------------------------------------------------------------- */

FixDivideTA::~FixDivideTA()
{
  delete random;
  delete [] group_id;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */

void FixDivideTA::post_constructor()
{
  // create fix nufeb/property/cycletime
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"div_ta_ct";
  fixarg[1] = group_id;
  fixarg[2] = (char *)"nufeb/property/cycletime";

  modify->add_fix(3, fixarg, 1);
  delete [] fixarg;
  fix_ct = (FixPropertyCycletime *)modify->fix[modify->nfix-1];
}

/* ---------------------------------------------------------------------- */

void FixDivideTA::init()
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double div_time = log(2) / 0.00000321 ;
      double prob = random->uniform();
      fix_ct->aprop[i][0] = prob * div_time;
      fix_ct->aprop[i][1] = prob * 216000;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDivideTA::compute()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;
  double **conc = grid->conc;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell =grid->cell(x[i]);
      double growth = grid->growth[igroup][cell][0];
      double div_time = log(2) / growth;

      // trigger cell division if current cell cycle time
      // is greater than calculated doubling time
      if (fix_ct->aprop[i][0] > div_time) {

        // attributes of the two daughter cells are identical to parent cell
        double theta = random->uniform() * 2 * MY_PI;
        double phi = random->uniform() * (MY_PI);

        double oldx = atom->x[i][0];
        double oldy = atom->x[i][1];
        double oldz = atom->x[i][2];

        // update daughter cell i
        double newx = oldx + (atom->shape[i][0] * cos(theta) * sin(phi) * DELTA);
        double newy = oldy + (atom->shape[i][1] * sin(theta) * sin(phi) * DELTA);
        //double newz = oldz + (atom->shape[i][2] * cos(phi) * DELTA*0.2);

        atom->x[i][0] = newx;
        atom->x[i][1] = newy;

        // create daughter cell j

        double *coord = new double[3];
        newx = oldx - (atom->shape[i][0] * cos(theta) * sin(phi) * DELTA);
        newy = oldy - (atom->shape[i][1] * sin(theta) * sin(phi) * DELTA);
        //newz = oldz + (atom->shape[i][2] * cos(phi) * DELTA*0.2);

        coord[0] = newx;
        coord[1] = newy;
        coord[2] = atom->x[i][2];

        atom->avec->create_atom(atom->type[i], coord);
        int j = atom->nlocal - 1;

        atom->tag[j] = 0;
        atom->mask[j] = atom->mask[i];
        atom->v[j][0] = atom->v[i][0];
        atom->v[j][1] = atom->v[i][1];
        atom->v[j][2] = atom->v[i][2];
        atom->f[j][0] = atom->f[i][0];
        atom->f[j][1] = atom->f[i][1];
        atom->f[j][2] = atom->f[i][2];

        atom->alpha[j] = atom->alpha[i];
        atom->shape[j][0] = atom->shape[i][0];
        atom->shape[j][1] = atom->shape[i][1];
        atom->shape[j][2] = atom->shape[i][2];
        atom->rmass[j] = atom->rmass[i];
        atom->biomass[j] = atom->biomass[i];


        modify->create_attribute(j);

        for (int m = 0; m < modify->nfix; m++)
          modify->fix[m]->update_arrays(i, j);

        delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}

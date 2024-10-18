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

#include "math_const.h"
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "compute.h"
#include "random_park.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "grid.h"
#include "group.h"
#include "fix_property_cycletime.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 1.05

/* ---------------------------------------------------------------------- */

FixDivideStem::FixDivideStem(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/division/stem requires skin atom style");

  if (narg < 6)
    error->all(FLERR, "Illegal fix epidermis/division/stem command");

  type_ta = utils::inumeric(FLERR,arg[3],true,lmp);
  mask_ta = group->find(arg[4]);

  if (mask_ta < 0)
    error->all(FLERR, "Can't find TA group in fix epidermis/division/stem command");
  mask_ta = 1 | group->bitmask[mask_ta];

  icyto = grid->find(arg[5]);
  if (icyto < 0)
    error->all(FLERR, "Can't find cytokine name");
  k_cyto = utils::numeric (FLERR,arg[6],true,lmp);
  if (k_cyto <= 0)
    error->all(FLERR, "cytokine affinity must be greater than zero");

  seed = utils::inumeric(FLERR,arg[7],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
  group_id = new char[strlen(arg[1])+1];
  strcpy(group_id, arg[1]);

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "pa") == 0) {
      pa1 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      pa2 = utils::numeric(FLERR,arg[iarg+2],true,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "pb") == 0 ) {
      pb1 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      pb2 = utils::numeric(FLERR,arg[iarg+2],true,lmp);
      iarg += 3;
    } else {
      error->all(FLERR, "Illegal fix epidermis/growth/diff command");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixDivideStem::~FixDivideStem()
{
  delete random;
  delete [] group_id;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */

void FixDivideStem::post_constructor()
{
  // create fix nufeb/property/cycletime
  char **fixarg = new char*[7];
  fixarg[0] = (char *)"div_sc_ct";
  fixarg[1] = group_id;
  fixarg[2] = (char *)"nufeb/property/cycletime";
  fixarg[3] = (char *)"seed";
  fixarg[4] = (char *)"2023";
  fixarg[5] = (char *)"max_time";
  fixarg[6] = (char *)"504000";

  modify->add_fix(7, fixarg, 1);
  delete [] fixarg;
  fix_ct = (FixPropertyCycletime *)modify->fix[modify->nfix-1];

}

/* ---------------------------------------------------------------------- */

void FixDivideStem::compute()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;
  double **conc = grid->conc;
  double *bulk = grid->bulk;

  int itype, jtype;
  int imask, jmask;
  int n = 0, flag = 1;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell =grid->cell(x[i]);
      double growth = grid->growth[igroup][cell][0];
      double div_time = log(2) / growth;
      //printf("SC %i = i %e %e\n", i, fix_ct->aprop[i][0], growth );
      // trigger cell division if current cell cycle time
      // is greater than calculated doubling time
      double ct = fix_ct->aprop[i][0];
      double pg = random->gaussian();

      if (ct + (ct*0.15*pg)  > div_time) {
        int cell = grid->cell(atom->x[i]);
        double prob = random->uniform();

        // apply hill function later
        double pa = pa1;
        double pb = pb1 - (pb2 * 1 / (1 + exp(k_cyto * (bulk[icyto] - 35))));
        if (flag == 1) {
          //if(comm->me == 0) printf("pa = %e pb = %e  f=%e\n", pa, pb, 1 / (1+exp(-k_cyto * (bulk[icyto] - 50))));
          flag = 0;
        }
        if (prob < pa) {
          // symmetric divide into two TA cells
          itype = type_ta;
          imask = mask_ta;
          jtype = type_ta;
          jmask = mask_ta;
        } else if (prob < 1-pb) {
          // self renewal
          itype = atom->type[i];
          imask = atom->mask[i];
          jtype = atom->type[i];
          jmask = atom->mask[i];
        } else {
          itype = atom->type[i];
          imask = atom->mask[i];
          jtype = type_ta;
          jmask = mask_ta;
        }

        double theta = random->uniform() * 2 * MY_PI;
        double phi = random->uniform() * (MY_PI);

        // attributes of the two daughter cells are identical to parent cell
        // update daughter cell i
        double oldx = atom->x[i][0];
        double oldy = atom->x[i][1];
        double oldz = atom->x[i][2];

        double newx = oldx + (atom->shape[i][0] * cos(theta) * sin(phi) * DELTA);
        double newy = oldy + (atom->shape[i][1] * sin(theta) * sin(phi) * DELTA);

        atom->x[i][0] = newx;
        atom->x[i][1] = newy;

        atom->type[i] = itype;
        atom->mask[i] = imask;

        // create daughter cell j

        double *coord = new double[3];
        newx = oldx - (atom->shape[i][0] * cos(theta) * sin(phi) * DELTA);
        newy = oldy - (atom->shape[i][1] * sin(theta) * sin(phi) * DELTA);

        coord[0] = newx;
        coord[1] = newy;
        coord[2] = atom->x[i][2];

        atom->avec->create_atom(jtype, coord);
        int j = atom->nlocal - 1;

        atom->tag[j] = 0;
        atom->mask[j] = jmask;
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
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_INT, MPI_SUM, world);

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

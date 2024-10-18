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
#include "update.h"
#include "error.h"
#include "domain.h"
#include "comm.h"

#include "fix_secretion.h"
#include "grid.h"
#include "group.h"
#include "modify.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSecretion::FixSecretion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR, "Illegal fix epidermis/growth/basal command");

  isub = -1;
  icyto = -1;
  cyto_max = 0.0;

  tstem = utils::inumeric (FLERR,arg[3],true,lmp);
  tta = utils::inumeric (FLERR,arg[4],true,lmp);

  isub = grid->find(arg[5]);
  if (isub < 0)
    error->all(FLERR, "Can't find growth factor name");

  icyto = grid->find(arg[6]);
  if (isub < 0)
    error->all(FLERR, "Can't find cytokine name");
  cyto_max = utils::numeric (FLERR,arg[7],true,lmp);
  if (cyto_max <= 0)
    error->all(FLERR, "cytokine affinity must be greater than zero");

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "yield") == 0 ) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0 ) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cyto_flag") == 0 ) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix epidermis/growth/stem command");
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixSecretion::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSecretion::biology_nufeb()
{
  secrete_gf();
  secrete_cyto();
}
/* ---------------------------------------------------------------------- */

void FixSecretion::secrete_gf()
{

  double **conc = grid->conc;
  double *bulk = grid->bulk;
  int nlocal = atom->nlocal;
  int nbasal = 0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->type[i] == tstem || atom->type[i] == tta) {
      nbasal++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &nbasal, 1, MPI_INT, MPI_SUM, world);

  bulk[isub] += (yield * nbasal - decay * bulk[isub]) * update->dt;

  //if(comm->me==0) printf("bulk-gf =%e\n", bulk[isub]);
}

/* ---------------------------------------------------------------------- */

void FixSecretion::secrete_cyto()
{
  double *bulk = grid->bulk;

  bulk[icyto] = cyto_max / (1 + exp(-0.01 * (bulk[isub] - 2500)));
  //if(comm->me==0) printf("bulk-cyto =%e\n", bulk[icyto]);
}
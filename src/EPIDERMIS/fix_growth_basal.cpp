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

#include "fix_growth_basal.h"
#include "grid.h"
#include "group.h"
#include "modify.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGrowthBasal::FixGrowthBasal(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR, "Illegal fix epidermis/growth/basal command");

  isub = -1;
  icyto = -1;
  growth = 0.0;
  yield = 1.0;
  decay = 0.0;
  cyto_affinity = 0.0;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find growth factor name");

  icyto = grid->find(arg[4]);
  if (icyto < 0)
    error->all(FLERR, "Can't find cytokine name");
  cyto_affinity = utils::numeric (FLERR,arg[5],true,lmp);
  if (cyto_affinity <= 0)
    error->all(FLERR, "cytokine affinity must be greater than zero");

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix epidermis/growth/stem command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthBasal::update_atoms()
{
  double *bulk = grid->bulk;
  int n;
  double m=0;
  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = 1;
    double tmp2 = 3 * bulk[icyto] / (cyto_affinity + bulk[icyto]);
//    n++;
//    m +=tmp2;
    grid->growth[igroup][i][0] = growth * (tmp1 + tmp2);
  }
  // printf("growth basal tmp2=%e \n", m/n);
}

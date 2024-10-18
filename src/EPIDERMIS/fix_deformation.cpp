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

#include "fix_deformation.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "grid.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "fix_property_diffstate.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeformation::FixDeformation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->skin_flag)
    error->all(FLERR, "fix epidermis/deformation requires skin atom style");

  if (narg < 4)
    error->all(FLERR,  "Illegal fix epidermis/deformation command");

  alpha_f = utils::numeric(FLERR,arg[3],true,lmp);
  p_f = utils::numeric(FLERR,arg[4],true,lmp);

  group_flag = 0;
  ngroups = 0;

  diffstate = nullptr;
  diffgroup = nullptr;

  memory->grow(diffstate,atom->ntypes+1,"fix skin/deformation:diffstate");
  memory->grow(diffgroup,atom->ntypes+1,"fix skin/deformation:diffgroup");

  for (int i = 1; i <= atom->ntypes; i++) {
    diffstate[i] = MAXFLOAT;
  }

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "adiff") == 0) {
      int igroup = group->find(arg[iarg+1]);
      int type = utils::inumeric(FLERR,arg[iarg+2],true,lmp);
      double state = utils::numeric(FLERR,arg[iarg+3],true,lmp);

      if (igroup == -1)
        error->all(FLERR, "Can't find group in fix skin/deformation command");

      diffstate[type] = state;
      diffgroup[type] = 1 | group->bitmask[igroup];
      iarg += 4;
    } else if (strcmp(arg[iarg], "groups") == 0) {
      group_flag = 1;
      ngroups = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      memory->create(groupbits,ngroups,"fix skin/deformation:groupbits");

      for (int i = 0; i < ngroups; i++) {
        int igroup = group->find(arg[iarg+2+i]);
        if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
        groupbits[i] = group->bitmask[igroup];
      }

      iarg = iarg + 2 + ngroups;
    } else {
      error->all(FLERR, "Illegal fix skin/deformation command");
    }
  }

  auto fixlist = modify->get_fix_by_style("^skin/property/diffstate");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/diffstate defined for this fix ");
  fix_diffstate = dynamic_cast<FixPropertyDiffstate *>(fixlist.front());
}


/* ---------------------------------------------------------------------- */

FixDeformation::~FixDeformation()
{
  memory->destroy(diffstate);
  memory->destroy(diffgroup);
}

/* ---------------------------------------------------------------------- */

int FixDeformation::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeformation::biology_nufeb()
{
  double ave1 = 0.0;
  double ave2 = 0.0;
  double ave3 = 0.0;
  double ave4 = 0.0;
  double ave5 = 0.0;
  double ave6 = 0.0;
  int n1 = 1;
  int n2 = 1;
  int n3 = 1;
  for (int i = 0; i < atom->nlocal; i++) {
    if (group_flag) {
      for (int j = 0; j < ngroups; j++) {
        if (atom->mask[i] & groupbits[j]) {
          compute(i);
          // printf("%i %e \n", i, fix_diffstate->vprop[i]);
        }
      }
    } else {
      if (atom->mask[i] & groupbit) {
        compute(i);
      }
    }
  }

//  for (int i = 0; i < atom->nlocal; i++) {
//    if (atom->type[i] == 4) {
//      ave1 += atom->shape[i][0]*2;
//      ave4 += atom->shape[i][2]*2;
//      n1++;
//    }
//    if (atom->type[i] == 5) {
//      ave2 += atom->shape[i][0]*2;
//      ave5 += atom->shape[i][2]*2;
//      n2++;
//    }
//    if (atom->type[i] == 6) {
//      ave3 += atom->shape[i][0]*2;
//      ave6 += atom->shape[i][2]*2;
//      n3++;
//    }
//  }
//  if (comm->me == 1)printf("spin ave = %e h=%e\n", ave1/n1, ave4/n1);
//  if (comm->me == 1) printf("gran ave = %e h=%e\n", ave2/n2, ave5/n2);
//  if (comm->me == 1) printf("corn ave = %e h = %e\n", ave3/n3, ave6/n3);
}

/* ---------------------------------------------------------------------- */

void FixDeformation::compute(int i)
{
  int ntypes = atom->ntypes;
  int *mask = atom->mask;

  double *alpha = atom->alpha;
  double **shape = atom->shape;
  int *type = atom->type;

  double tmp = std::pow(fix_diffstate->vprop[i]/3.0, p_f);
  alpha[i] = (1 + alpha_f *  tmp) / (1 + tmp);
  shape[i][0] = shape[i][0] * alpha[i];
  shape[i][1] = shape[i][1] * alpha[i];
  shape[i][2] = shape[i][2] / (alpha[i] * alpha[i]);

  for (int j = 1; j <= ntypes; j++) {
    if (fix_diffstate->vprop[i] > diffstate[j]) {
      type[i] = j;
      mask[i] = diffgroup[j];
    }
  }
}

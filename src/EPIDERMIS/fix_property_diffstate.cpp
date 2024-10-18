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

#include "fix_property_diffstate.h"

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

FixPropertyDiffstate::FixPropertyDiffstate(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/diffstate command");

  create_attribute = 1;
  // use vprop if size_peratom_cols = 1
  size_peratom_cols = 1;
  grow_arrays(atom->nmax);

  ical = -1;
  ks = 0.0;
  group_flag = 0;
  ngroups = 0;

  ical = grid->find(arg[3]);
  if (ical < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "r_max") == 0) {
      r_max = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if(strcmp(arg[iarg], "ks") == 0) {
      ks = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "groups") == 0) {
      group_flag = 1;
      ngroups = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      memory->create(groupbits,ngroups,"fix epidermis/differentiation:groupbits");

      for (int i = 0; i < ngroups; i++) {
        int igroup = group->find(arg[iarg+2+i]);
        if (igroup == -1) error->all(FLERR,"Could not find fix group ID");
        groupbits[i] = group->bitmask[igroup];
      }

      iarg = iarg + 2 + ngroups;
    } else {
      error->all(FLERR, "Illegal fix nufeb/property/diffstate command");
    }
  }

  for (int i = 0; i < atom->nlocal; i++) {
    vprop[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */
FixPropertyDiffstate::~FixPropertyDiffstate()
{
  if (group_flag)
    memory->destroy(groupbits);
}

/* ---------------------------------------------------------------------- */

int FixPropertyDiffstate::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropertyDiffstate::biology_nufeb()
{
  double **conc = grid->conc;

//  std::vector<double> v;
//  double total = 0.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (group_flag) {
      for (int j = 0; j < ngroups; j++){
        if (atom->mask[i] & groupbits[j]) {
          int cell = grid->cell(atom->x[i]);
          double r = r_max / (1 + exp(-500 * (conc[ical][cell] - ks)));
          vprop[i] = vprop[i] + r * update->dt;
//          if (conc[ical][cell] > 0.03 && conc[ical][cell] < 0.033) {
//            v.push_back(vprop[i]);
//            total += vprop[i];
//          }
          //if (conc[ical][cell] > 0.03 && conc[ical][cell] < 0.033) printf("r = %e v=%e \n",r, vprop[i]);
        }
      }
    } else {
      if (atom->mask[i] & groupbit) {
        int cell = grid->cell(atom->x[i]);
        double r = r_max / (1 + exp(-500 * (conc[ical][cell] - ks)));
        vprop[i] = vprop[i] + r * update->dt;
      }
    }
  }
//  double mean = total/v.size();
//  double variance = 0.0;
//
//  for (int i = 0; i < v.size(); i++) {
//    variance += pow(v[i] - mean, 2);
//  }
//
//  variance /= v.size();
  //printf("ave = %e std = %e \n",mean, sqrt(variance));
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyDiffstate::set_arrays(int j)
{
  vprop[j] = 0.0;
}

/* ----------------------------------------------------------------------
   update array values of two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyDiffstate::update_arrays(int i, int j)
{
  vprop[j] = 0.0;
  vprop[i] = 0.0;
}



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

#include "compute_nufeb_biomass.h"

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebBiomass::ComputeNufebBiomass(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute biomass command");

  vector_flag = 1;
  size_vector = atom->ntypes+1;
  extvector = 0;

  vector = new double[atom->ntypes+1]();
}

/* ---------------------------------------------------------------------- */

ComputeNufebBiomass::~ComputeNufebBiomass()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeNufebBiomass::compute_vector()
{
  invoked_vector = update->ntimestep;

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  for (int i = 0; i < ntypes + 1; i++) {
    vector[i] = 0;
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      int t = type[i] - 1;
      vector[t] += rmass[i];
      vector[ntypes] += rmass[i];
    }
}
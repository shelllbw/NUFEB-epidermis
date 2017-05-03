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

#include "fix_bio_death.h"

#include <string.h>

#include "atom.h"
#include "atom_vec.h"
#include "error.h"

#include "atom_vec_bio.h"
#include "force.h"
#include "input.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeath::FixDeath(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix death requires atom style bio");

  if (narg != 5) error->all(FLERR,"Illegal fix death command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix death command");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[4][2]);

  //force_reneighbor = 1;
  //next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixDeath::~FixDeath(){
  delete [] var;
}

/* ---------------------------------------------------------------------- */

int FixDeath::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeath::init(){
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix death does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix death is invalid style");

  deadMass= input->variable->compute_equal(ivar);

  if (avec->typeDEAD == 0) {
    error->all(FLERR,"At least one initial DEAD bactera is required.");
  }
}

/* ---------------------------------------------------------------------- */

void FixDeath::pre_exchange()
{
  //if (next_reneighbor != update->ntimestep) return;
  death();
}

/* ---------------------------------------------------------------------- */

void FixDeath::death()
{
  //double criticalMass = 1e-20;
  int *type = atom->type;
  int *mask = atom->mask;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

//(rmass[i] < criticalMass)
  for (i = nall-1; i >= 0; i--) {
//  	//delete atom
//  	if((mask[i] & groupbit) && (rmass[i] < criticalMass)) {
//			atom->avec->copy(nall-1,i,1);
//			atom->nlocal--;
//			atom->natoms--;
//			continue;
//  	}

    if ((mask[i] & groupbit) && (mask[i] != avec->maskDEAD) && (mask[i] != avec->maskEPS)) {
    	if(rmass[i] < deadMass) {
    		type[i] = avec->typeDEAD;
    		mask[i] = avec->maskDEAD;
    	}
    }
  }

//  if (atom->map_style) {
//    atom->nghost = 0;
//    atom->map_init();
//    atom->map_set();
//  }
 // printf("nlocal=%i, all=%i \n", nlocal, atom->natoms );
  //next_reneighbor += nevery;
}
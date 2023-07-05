/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "atom_vec_skin.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecSkin::AtomVecSkin(LAMMPS *lmp) : AtomVec(lmp)
{
  mass_type = PER_ATOM;
  molecular = Atom::ATOMIC;

  atom->skin_flag = 1;
  atom->rmass_flag = 1;
  atom->biomass_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"shape", "rmass", "biomass"};
  fields_copy = {"shape", "rmass", "biomass"};
  fields_comm = {""};
  fields_comm_vel = {""};
  fields_reverse = {""};
  fields_border = {"shape", "rmass", "biomass"};
  fields_border_vel = {"shape", "rmass", "biomass"};
  fields_exchange = {"shape", "rmass", "biomass"};
  fields_restart = {"shape", "rmass", "biomass"};
  fields_create = {"shape", "rmass", "biomass"};
  fields_data_atom = {"id", "type", "shape", "rmass", "x"};
  fields_data_vel =  {"id", "v"};
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecSkin::process_args(int narg, char **arg)
{
  if (narg != 0 && narg != 1)
    error->all(FLERR,"Illegal atom_style skin command");

  shapevary = 1;
  if (narg == 1) {
    shapevary = utils::numeric(FLERR,arg[0],true,lmp);
    if (shapevary < 0 || shapevary > 1)
      error->all(FLERR,"Illegal atom_style skin command");
  }

  // dynamic particle properties must be communicated every step

  if (shapevary) {
    fields_comm = {"shape", "rmass"};
    fields_comm_vel = {"shape", "rmass"};
  }

  // delay setting up of fields until now

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecSkin::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag && shapevary == 0)
        error->all(FLERR,"Fix adapt changes particle radii "
                         "but atom_style coccus is not dynamic");
    } else if (strcmp(modify->fix[i]->style,"nufeb/monod") == 0) {
      if (shapevary == 0)
        error->all(FLERR,"Fix nufeb/monod changes particle radii "
                         "but atom_style coccus is not dynamic");
    } else if (strcmp(modify->fix[i]->style,"nufeb/divide") == 0) {
      if (shapevary == 0)
        error->all(FLERR,"Fix nufeb/divide changes particle radii "
                         "but atom_style coccus is not dynamic");
    }
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSkin::grow_pointers()
{
  shape = atom->shape;
  rmass = atom->rmass;
  biomass = atom->biomass;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecSkin::create_atom_post(int ilocal)
{
  shape[ilocal][0] = shape[ilocal][1] = shape[ilocal][2] = 1e-5;
  rmass[ilocal] = 4.0*MY_PI/3.0 * 1e-5*1e-5*1e-5;
  biomass[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSkin::data_atom_post(int ilocal)
{
  radius_one = 0.5 * atom->radius[ilocal];
  radius[ilocal] = radius_one;
  if (radius_one > 0.0)
    rmass[ilocal] *= 4.0*MY_PI/3.0 * radius_one*radius_one*radius_one;

  if (rmass[ilocal] <= 0.0)
    error->one(FLERR,"Invalid density in Atoms section of data file");

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;

  outer_radius_one = 0.5 * atom->outer_radius[ilocal];
  outer_radius[ilocal] = outer_radius_one;
  if (outer_radius[ilocal] < radius[ilocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }

  outer_mass[ilocal] = (4.0*MY_PI/3.0)*
                       ((outer_radius[ilocal]*outer_radius[ilocal]*outer_radius[ilocal])
                        -(radius[ilocal]*radius[ilocal]*radius[ilocal])) * 30;

  biomass[ilocal] = 1.0;
}
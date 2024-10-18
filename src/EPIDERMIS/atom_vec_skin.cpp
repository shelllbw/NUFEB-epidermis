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

  fields_grow = {"shape", "rmass", "biomass", "alpha"};
  fields_copy = {"shape", "rmass", "biomass", "alpha"};
  fields_border = {"shape", "rmass", "biomass", "alpha"};
  fields_border_vel = {"shape", "rmass", "biomass", "alpha"};
  fields_exchange = {"shape", "rmass", "biomass", "alpha"};
  fields_restart = {"shape", "rmass", "biomass", "alpha"};
  fields_create = {"shape", "rmass", "biomass", "alpha"};
  fields_data_atom = {"id", "type", "rmass", "x", "shape"};
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

  // check if optional shapevary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag && shapevary == 0)
        error->all(FLERR,"Fix adapt changes particle radii "
                         "but atom_style skin is not dynamic");
    } else if (strcmp(modify->fix[i]->style,"nufeb/division") == 0) {
      if (shapevary == 0)
        error->all(FLERR,"Fix nufeb/divide changes particle radii "
                         "but atom_style skin is not dynamic");
    }
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSkin::grow_pointers()
{
  shape = atom->shape;
  alpha = atom->alpha;
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
  alpha[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSkin::data_atom_post(int ilocal)
{
  shape_a = 0.5 * shape[ilocal][0];
  shape_b = 0.5 * shape[ilocal][1];
  shape_c = 0.5 * shape[ilocal][2];

  shape[ilocal][0] = shape_a;
  shape[ilocal][1] = shape_b;
  shape[ilocal][2] = shape_c;

  if (shape_a<=0 || shape_b<=0 || shape_c<=0)
    error->one(FLERR,"Invalid cell shape in Atoms section of data file");

  rmass[ilocal] *= 4.0*MY_PI/3.0 * shape_a*shape_b*shape_c;

  biomass[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecSkin::pack_data_pre(int ilocal)
{
  shape_a = shape[ilocal][0];
  shape_b = shape[ilocal][1];
  shape_c = shape[ilocal][2];
  rmass_one = rmass[ilocal];
  alpha_one = alpha[ilocal];

  shape[ilocal][0] *= 2.0;
  shape[ilocal][1] *= 2.0;
  shape[ilocal][2] *= 2.0;

  if (shape_a > 0 && shape_b >0 && shape_c >0)
    rmass[ilocal] =
        rmass_one / (4.0*MY_PI/3.0 * shape_a*shape_b*shape_c);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecSkin::pack_data_post(int ilocal)
{
  shape[ilocal][0] = shape_a;
  shape[ilocal][1] = shape_b;
  shape[ilocal][2] = shape_c;
  rmass[ilocal] = rmass_one;
}

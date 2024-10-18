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
#include <cmath>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "variable.h"
#include "comm.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "math_const.h"
#include "group.h"

#include "fix_update_base.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static constexpr double BIG = 1.0e30;
static constexpr double EPSILON = 1.0e-6;
static constexpr double LB_FACTOR = 1.1;

enum { COUNT, INSERT, INSERT_SELECTED };

/* ---------------------------------------------------------------------- */

FixUpdateBase::FixUpdateBase(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // check for compatible lattice
  int latsty = domain->lattice->style;
  if (latsty == Lattice::SQ || latsty == Lattice::SQ2 || latsty == Lattice::HEX)
    error->all(FLERR, "Lattice style incompatible with simulation dimension");

  varflag = 0;
  vstr = xstr = ystr = zstr = nullptr;
  dlist = nullptr;

  ntype = utils::inumeric(FLERR, arg[3], false, lmp);
  if ((ntype <= 0) || (ntype > atom->ntypes))
    error->all(FLERR, "Invalid atom type in create_atoms command");

  mask_base = 1 | group->bitmask[igroup];
  nbasis = domain->lattice->nbasis;
  basistype = new int[nbasis];
  for (int i = 0; i < nbasis; i++) basistype[i] = ntype;

  int iarg;
  if (strcmp(arg[4], "region") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "fix skin/update_base region", error);
    region = domain->get_region_by_id(arg[5]);
    if (!region) error->all(FLERR, "fix skin/update_base region {} does not exist", arg[2]);
    region->init();
    region->prematch();
    iarg = 6;
  } else
    error->all(FLERR, "Unknown fix skin/update_base command option {}", arg[1]);

  while (iarg < narg) {
    if (strcmp(arg[iarg], "var") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "create_atoms var", error);
      delete[] vstr;
      vstr = utils::strdup(arg[iarg + 1]);
      varflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "set") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "create_atoms set", error);
      if (strcmp(arg[iarg + 1], "x") == 0) {
        delete[] xstr;
        xstr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "y") == 0) {
        delete[] ystr;
        ystr = utils::strdup(arg[iarg + 2]);
      } else if (strcmp(arg[iarg + 1], "z") == 0) {
        delete[] zstr;
        zstr = utils::strdup(arg[iarg + 2]);
      } else
        error->all(FLERR, "Unknown fix skin/update_base set option {}", arg[iarg + 2]);
      iarg += 3;
    }
  }

  if (!vstr && (xstr || ystr || zstr))
    error->all(FLERR, "Incomplete use of variables in  fix skin/update_base command");
  if (vstr && (!xstr && !ystr && !zstr))
    error->all(FLERR, "Incomplete use of variables in  fix skin/update_base command");

  if (varflag) {
    vvar = input->variable->find(vstr);

    if (vvar < 0) error->all(FLERR, "Variable {} for fix skin/update_base does not exist", vstr);
    if (!input->variable->equalstyle(vvar))
      error->all(FLERR, "Variable for fix skin/update_base is invalid style");

    if (xstr) {
      xvar = input->variable->find(xstr);
      if (xvar < 0) error->all(FLERR, "Variable {} for fix skin/update_base does not exist", xstr);
      if (!input->variable->internalstyle(xvar))
        error->all(FLERR, "Variable for fix skin/update_base is invalid style");
    }
    if (ystr) {
      yvar = input->variable->find(ystr);
      if (yvar < 0) error->all(FLERR, "Variable {} for fix skin/update_base does not exist", ystr);
      if (!input->variable->internalstyle(yvar))
        error->all(FLERR, "Variable for fix skin/update_base is invalid style");
    }
    if (zstr) {
      zvar = input->variable->find(zstr);
      if (zvar < 0) error->all(FLERR, "Variable {} for fix skin/update_base does not exist", zstr);
      if (!input->variable->internalstyle(zvar))
        error->all(FLERR, "Variable for fix skin/update_base is invalid style");
    }
  }

  if (nbasis == 0) error->all(FLERR, "Cannot create atoms with undefined lattice");

  triclinic = domain->triclinic;

  double epsilon[3];
  if (triclinic)
    epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
  else {
    epsilon[0] = domain->prd[0] * EPSILON;
    epsilon[1] = domain->prd[1] * EPSILON;
    epsilon[2] = domain->prd[2] * EPSILON;
  }

  if (triclinic == 0) {
    sublo[0] = domain->sublo[0];
    subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1];
    subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2];
    subhi[2] = domain->subhi[2];
  } else {
    sublo[0] = domain->sublo_lamda[0];
    subhi[0] = domain->subhi_lamda[0];
    sublo[1] = domain->sublo_lamda[1];
    subhi[1] = domain->subhi_lamda[1];
    sublo[2] = domain->sublo_lamda[2];
    subhi[2] = domain->subhi_lamda[2];
  }

  if (comm->layout != Comm::LAYOUT_TILED) {
    if (domain->xperiodic) {
      if (comm->myloc[0] == 0) sublo[0] -= epsilon[0];
      if (comm->myloc[0] == comm->procgrid[0] - 1) subhi[0] -= 2.0 * epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->myloc[1] == 0) sublo[1] -= epsilon[1];
      if (comm->myloc[1] == comm->procgrid[1] - 1) subhi[1] -= 2.0 * epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->myloc[2] == 0) sublo[2] -= epsilon[2];
      if (comm->myloc[2] == comm->procgrid[2] - 1) subhi[2] -= 2.0 * epsilon[2];
    }
  } else {
    if (domain->xperiodic) {
      if (comm->mysplit[0][0] == 0.0) sublo[0] -= epsilon[0];
      if (comm->mysplit[0][1] == 1.0) subhi[0] -= 2.0 * epsilon[0];
    }
    if (domain->yperiodic) {
      if (comm->mysplit[1][0] == 0.0) sublo[1] -= epsilon[1];
      if (comm->mysplit[1][1] == 1.0) subhi[1] -= 2.0 * epsilon[1];
    }
    if (domain->zperiodic) {
      if (comm->mysplit[2][0] == 0.0) sublo[2] -= epsilon[2];
      if (comm->mysplit[2][1] == 1.0) subhi[2] -= 2.0 * epsilon[2];
    }
  }
}

/* ---------------------------------------------------------------------- */
FixUpdateBase::~FixUpdateBase()
{
  delete[] basistype;
  delete[] vstr;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
}

/* ---------------------------------------------------------------------- */

int FixUpdateBase::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUpdateBase::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  // record wall time for atom creation
  delete_base();

  MPI_Barrier(world);
  double time1 = platform::walltime();

  // clear ghost count and any ghost bonus data internal to AtomVec
  // same logic as beginning of Comm::exchange()
  // do it now b/c creating atoms will overwrite ghost atoms

  atom->nghost = 0;
  atom->avec->clear_bonus();

  // add atoms/molecules in one of 3 ways

  bigint natoms_previous = atom->natoms;
  int nlocal_previous = atom->nlocal;

  add_lattice();

  // init per-atom fix/compute/variable values for created atoms

  atom->data_fix_compute_variable(nlocal_previous, atom->nlocal);

  // set new total # of atoms and error check

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT) error->all(FLERR, "Too many total atoms");

  // add IDs for newly created atoms
  // check that atom IDs are valid

  if (atom->tag_enable) atom->tag_extend();
  atom->tag_check();

  // if global map exists, reset it
  // invoke map_init() b/c atom count has grown

  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }

  // print status

//  MPI_Barrier(world);
//  if (comm->me == 0) {
//    utils::logmesg(lmp, "Created {} atoms\n", atom->natoms - natoms_previous);
//
//    domain->print_box("  using box units in ");
//    utils::logmesg(lmp, "  create_atoms CPU = {:.3f} seconds\n", platform::walltime() - time1);
//  }
}


/* ---------------------------------------------------------------------- */

void FixUpdateBase::add_lattice()
{
  // convert 8 corners of my subdomain from box coords to lattice coords
  // for orthogonal, use corner pts of my subbox
  // for triclinic, use bounding box of my subbox
  // xyz min to max = bounding box around the domain corners in lattice space

  double bboxlo[3], bboxhi[3];

  if (triclinic == 0) {
    bboxlo[0] = domain->sublo[0];
    bboxhi[0] = domain->subhi[0];
    bboxlo[1] = domain->sublo[1];
    bboxhi[1] = domain->subhi[1];
    bboxlo[2] = domain->sublo[2];
    bboxhi[2] = domain->subhi[2];
  } else
    domain->bbox(domain->sublo_lamda, domain->subhi_lamda, bboxlo, bboxhi);

  // narrow down the subbox by the bounding box of the given region, if available.
  // for small regions in large boxes, this can result in a significant speedup

  if (region->bboxflag) {

    const double rxmin = region->extent_xlo;
    const double rxmax = region->extent_xhi;
    const double rymin = region->extent_ylo;
    const double rymax = region->extent_yhi;
    const double rzmin = region->extent_zlo;
    const double rzmax = region->extent_zhi;

    if (rxmin > bboxlo[0]) bboxlo[0] = (rxmin > bboxhi[0]) ? bboxhi[0] : rxmin;
    if (rxmax < bboxhi[0]) bboxhi[0] = (rxmax < bboxlo[0]) ? bboxlo[0] : rxmax;
    if (rymin > bboxlo[1]) bboxlo[1] = (rymin > bboxhi[1]) ? bboxhi[1] : rymin;
    if (rymax < bboxhi[1]) bboxhi[1] = (rymax < bboxlo[1]) ? bboxlo[1] : rymax;
    if (rzmin > bboxlo[2]) bboxlo[2] = (rzmin > bboxhi[2]) ? bboxhi[2] : rzmin;
    if (rzmax < bboxhi[2]) bboxhi[2] = (rzmax < bboxlo[2]) ? bboxlo[2] : rzmax;
  }

  double xmin, ymin, zmin, xmax, ymax, zmax;
  xmin = ymin = zmin = BIG;
  xmax = ymax = zmax = -BIG;

  // convert to lattice coordinates and set bounding box

  domain->lattice->bbox(1, bboxlo[0], bboxlo[1], bboxlo[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxhi[0], bboxlo[1], bboxlo[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxlo[0], bboxhi[1], bboxlo[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxhi[0], bboxhi[1], bboxlo[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxlo[0], bboxlo[1], bboxhi[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxhi[0], bboxlo[1], bboxhi[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxlo[0], bboxhi[1], bboxhi[2], xmin, ymin, zmin, xmax, ymax, zmax);
  domain->lattice->bbox(1, bboxhi[0], bboxhi[1], bboxhi[2], xmin, ymin, zmin, xmax, ymax, zmax);

  // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my subbox
  // overlap = any part of a unit cell (face,edge,pt) in common with my subbox
  // in lattice space, subbox is a tilted box
  // but bbox of subbox is aligned with lattice axes
  // so ilo:khi unit cells should completely tile bounding box
  // decrement lo, increment hi to avoid round-off issues in lattice->bbox(),
  //   which can lead to missing atoms in rare cases
  // extra decrement of lo if min < 0, since static_cast(-1.5) = -1

  ilo = static_cast<int>(xmin) - 1;
  jlo = static_cast<int>(ymin) - 1;
  klo = static_cast<int>(zmin) - 1;
  ihi = static_cast<int>(xmax) + 1;
  jhi = static_cast<int>(ymax) + 1;
  khi = static_cast<int>(zmax) + 1;

  if (xmin < 0.0) ilo--;
  if (ymin < 0.0) jlo--;
  if (zmin < 0.0) klo--;

  // count lattice sites on each proc

  nlatt_overflow = 0;
  loop_lattice(COUNT);

  // nadd = # of atoms each proc will insert (estimated if subsetflag)

  int overflow;
  MPI_Allreduce(&nlatt_overflow, &overflow, 1, MPI_INT, MPI_SUM, world);
  if (overflow) error->all(FLERR, "Create_atoms lattice size overflow on 1 or more procs");

  bigint nadd;

  if (comm->nprocs == 1)
    nadd = nlatt;
  else
    nadd = static_cast<bigint>(LB_FACTOR * nlatt);

  // allocate atom arrays to size N, rounded up by AtomVec->DELTA

  bigint nbig = atom->avec->roundup(nadd + atom->nlocal);
  int n = static_cast<int>(nbig);
  atom->avec->grow(n);

  // add atoms or molecules
  // if no subset: add to all lattice sites
  // if subset: count lattice sites, select random subset, then add

  loop_lattice(INSERT);

}

/* ----------------------------------------------------------------------
   iterate on 3d periodic lattice of unit cells using loop bounds
   iterate on nbasis atoms in each unit cell
   convert lattice coords to box coords
   check if lattice point meets all criteria to be added
   perform action on atom or molecule (on each basis point) if meets all criteria
   actions = add, count, add if flagged
------------------------------------------------------------------------- */

void FixUpdateBase::loop_lattice(int action)
{
  int i, j, k, m;

  const double *const *const basis = domain->lattice->basis;

  nlatt = 0;

  for (k = klo; k <= khi; k++) {
    for (j = jlo; j <= jhi; j++) {
      for (i = ilo; i <= ihi; i++) {
        for (m = 0; m < nbasis; m++) {
          double *coord;
          double x[3], lamda[3];

          x[0] = i + basis[m][0];
          x[1] = j + basis[m][1];
          x[2] = k + basis[m][2];

          // convert from lattice coords to box coords

          domain->lattice->lattice2box(x[0], x[1], x[2]);

          // if a region was specified, test if atom is in it
          if (!region->match(x[0], x[1], x[2])) continue;

          // if variable test specified, eval variable

          if (varflag && vartest(x) == 0) continue;

          // test if atom/molecule position is in my subbox

          if (triclinic) {
            domain->x2lamda(x, lamda);
            coord = lamda;
          } else
            coord = x;

          if (coord[0] < sublo[0] || coord[0] >= subhi[0] || coord[1] < sublo[1] ||
              coord[1] >= subhi[1] || coord[2] < sublo[2] || coord[2] >= subhi[2])
            continue;

          // this proc owns the lattice site
          // perform action: add, just count, add if flagged
          // add = add an atom or entire molecule to my list of atoms

          if (action == INSERT) {
            atom->avec->create_atom(basistype[m], x);
            int n = atom->nlocal - 1;

            atom->tag[n] = 0;
            atom->mask[n] = mask_base;
            atom->v[n][0] = 0.0;
            atom->v[n][1] = 0.0;
            atom->v[n][2] = 0.0;
            atom->f[n][0] = 0.0;
            atom->f[n][1] = 0.0;
            atom->f[n][2] = 0.0;

            atom->alpha[n] = 0.0;
            atom->shape[n][0] = 5e-6;
            atom->shape[n][1] = 5e-6;
            atom->shape[n][2] = 5e-6;
            atom->rmass[i] = 4.0 * MY_PI / 3.0 *
                             atom->shape[i][0] * atom->shape[i][1] * atom->shape[i][0] * 4000;
            atom->biomass[n] = 1.0;

            modify->create_attribute(n);

          } else if (action == COUNT) {
            if (nlatt == MAXSMALLINT) nlatt_overflow = 1;
          }

          nlatt++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   test a generated atom position against variable evaluation
   first set x,y,z values in internal variables
------------------------------------------------------------------------- */

int FixUpdateBase::vartest(double *x)
{
  if (xstr) input->variable->internal_set(xvar, x[0]);
  if (ystr) input->variable->internal_set(yvar, x[1]);
  if (zstr) input->variable->internal_set(zvar, x[2]);

  double value = input->variable->compute_equal(vvar);

  if (value == 0.0) return 0;
  return 1;
}


void FixUpdateBase::delete_base()
{
  // allocate and initialize deletion list

  bigint natoms_previous = atom->natoms;

  int nlocal = atom->nlocal;
  memory->create(dlist, nlocal, "delete_atoms:dlist");
  for (int i = 0; i < nlocal; i++) dlist[i] = 0;

  int *mask = atom->mask;
  int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dlist[i] = 1;

  // delete local atoms flagged in dlist
  // reset nlocal

  AtomVec *avec = atom->avec;

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal - 1, i, 1);
      dlist[i] = dlist[nlocal - 1];
      nlocal--;
    } else
      i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  // reset atom->natoms and also topology counts

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);

  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped

  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

//  bigint ndelete = natoms_previous - atom->natoms;
//  if (comm->me == 0) {
//    printf("Deleted %i atoms \n", ndelete);
//  }
}
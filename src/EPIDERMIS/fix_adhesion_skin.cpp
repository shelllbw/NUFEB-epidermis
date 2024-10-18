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

#include "fix_adhesion_skin.h"

#include <math.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdhesionSkin::FixAdhesionSkin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix skin/adhesion command");

  virial_global_flag = 1;
  allocated = 0;

  af = nullptr;

  zeta = utils::numeric(FLERR,arg[3],true,lmp);
  r = utils::numeric(FLERR,arg[4],true,lmp);

  if (zeta <= 0)
    error->all(FLERR, "Illegal fix skin/adhesion command: zeta must be great than zero");
}

/* ---------------------------------------------------------------------- */

FixAdhesionSkin::~FixAdhesionSkin()
{
  memory->destroy(af);
}

/* ---------------------------------------------------------------------- */

int FixAdhesionSkin::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "af") == 0) {
      if (!allocated) allocate();

      int ilo,ihi,jlo,jhi;
      utils::bounds(FLERR,arg[iarg+1],1,atom->ntypes,ilo,ihi,error);
      utils::bounds(FLERR,arg[iarg+2],1,atom->ntypes,jlo,jhi,error);

      double sigma_one = utils::numeric(FLERR,arg[iarg+3],true,lmp);

      int count = 0;
      for (int i = ilo; i <= ihi; i++) {
        for (int j = MAX(jlo,i); j <= jhi; j++) {
          af[i][j] = sigma_one;
          af[j][i] = sigma_one;
          count++;
        }
      }
      if (count == 0) error->all(FLERR,"Incorrect args for af coefficient");
      iarg += 4;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void FixAdhesionSkin::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(af,atom->ntypes+1,atom->ntypes+1,"nufeb/cohesion:ah");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      af[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixAdhesionSkin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdhesionSkin::init() {
  if (!allocated) error->all(FLERR,"fix adhesion coeffs are not set");

  neighbor->add_request(this,NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

void FixAdhesionSkin::init_list(int id, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixAdhesionSkin::post_force(int vflag)
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixAdhesionSkin::compute()
{
  double minw, mind, minh;

  int nlocal = atom->nlocal;
  int *type = atom->type;
  double **x = atom->x;
  double **f = atom->f;
  double **shape = atom->shape;
  int newton_pair = force->newton_pair;
  
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int i = 0; i < 6; i++)
    virial[i] = 0.0;
  
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double xtmp = x[i][0];
    double ytmp = x[i][1];
    double ztmp = x[i][2];
    double xshape = shape[i][0];
    double yshape = shape[i][1];
    double zshape = shape[i][2];

    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      int jtype = type[j];

      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];

      double sumw = xshape + shape[j][0];
      double sumd = yshape + shape[j][1];
      double sumh = zshape + shape[j][2];

      minw = MAX(fabs(delx)-sumw, 0);
      mind = MAX(fabs(dely)-sumd, 0);
      minh = MAX(fabs(delz)-sumh, 0);

      double gij = sqrt(minw*minw + mind*mind + minh*minh);
      double rij = gij/r;

      double c1 = 1/ sqrt(2*zeta);
      double c2 = c1* exp(-0.5);

      double afij = af[itype][jtype];
      double ccel = -afij*((rij + c1
          )*exp(-zeta*(rij+c1)*(rij+c1)) - c2*exp(-zeta*rij*rij));

      double ccelx = delx*ccel;
      double ccely = dely*ccel;
      double ccelz = delz*ccel;

      f[i][0] -= ccelx;
      f[i][1] -= ccely;
      f[i][2] -= ccelz;

      if (newton_pair || j < nlocal){
        f[j][0] += ccelx;
        f[j][1] += ccely;
        f[j][2] += ccelz;
      }
    }
  }
}

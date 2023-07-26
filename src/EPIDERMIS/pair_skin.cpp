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

#include "pair_skin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSkin::PairSkin(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairSkin::~PairSkin ()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(k_n);
    memory->destroy(setflag);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSkin::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSkin::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSkin::coeff (int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double kn_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      k_n[i][j] = kn_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for ah coefficient");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSkin::init_one (int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  k_n[j][i] = k_n[i][j];
  return cut_global;
}


/* ---------------------------------------------------------------------- */

void PairSkin::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fx,fy,fz;
  double xshape, yshape, zshape;
  double sumx, sumy, sumz;
  double minw, mind, minh, sij;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **shape = atom->shape;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      xshape = shape[i][0];
      yshape = shape[i][1];
      zshape = shape[i][2];
      itype = type[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        sumx = (xshape + shape[j][0]) * 0.5;
        sumy = (yshape + shape[j][1]) * 0.5;
        sumz = (zshape + shape[j][2]) * 0.5;

        // overlapping area between two ellipsoids is
        // approximated by the overlap between the two
        // cells as if they are rectangles with widths and heights identical
        // to the corresponding ellipsoids.

        if (abs(delx) < sumx || abs(dely) < sumy || abs(delz) < sumz)
          {

            minw = MIN(sumx - delx, xshape);
            minw = MIN(minw, shape[j][0]);
            mind = MIN(sumy - dely, yshape);
            mind = MIN(mind, shape[j][1]);
            minh = MIN(sumz - delz, zshape);
            minh = MIN(minh, shape[j][2]);

            // overlap area
            sij = minw * mind * minh;

            fx = delx * sij * k_n[i][j];
            fy = dely * sij * k_n[i][j];
            fz = delz * sij * k_n[i][j];

            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;

            if (newton_pair || j < nlocal)
              {
                f[j][0] -= fx;
                f[j][1] -= fy;
                f[j][2] -= fz;
              }
          }
        }

      }
  }
}


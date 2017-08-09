/*
 * fix_kinetics/diffusionS.cpp
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "atom.h"
#include "domain.h"
#include "error.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "atom_vec_bio.h"
#include "fix_bio_kinetics_diffusionS.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "memory.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsDiffusionS::FixKineticsDiffusionS(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 11) error->all(FLERR,"Not enough arguments in fix diffusion command");

  var = new char*[4];
  ivar = new int[4];

  if(strcmp(arg[3], "exp") == 0) sflag = 0;
  else if(strcmp(arg[3], "imp") == 0) sflag = 1;
  else error->all(FLERR,"Illegal PDE method command");

  for (int i = 0; i < 4; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if(strcmp(arg[8], "pp") == 0) xbcflag = 0;
  else if(strcmp(arg[8], "dd") == 0) xbcflag = 1;
  else if(strcmp(arg[8], "nd") == 0) xbcflag = 2;
  else if(strcmp(arg[8], "nn") == 0) xbcflag = 3;
  else if(strcmp(arg[8], "dn") == 0) xbcflag = 4;
  else error->all(FLERR,"Illegal x-axis boundary condition command");

  if(strcmp(arg[9], "pp") == 0) ybcflag = 0;
  else if(strcmp(arg[9], "dd") == 0) ybcflag = 1;
  else if(strcmp(arg[9], "nd") == 0) ybcflag = 2;
  else if(strcmp(arg[9], "nn") == 0) ybcflag = 3;
  else if(strcmp(arg[9], "dn") == 0) ybcflag = 4;
  else error->all(FLERR,"Illegal y-axis boundary condition command");

  if(strcmp(arg[10], "pp") == 0) zbcflag = 0;
  else if(strcmp(arg[10], "dd") == 0) zbcflag = 1;
  else if(strcmp(arg[10], "nd") == 0) zbcflag = 2;
  else if(strcmp(arg[10], "nn") == 0) zbcflag = 3;
  else if(strcmp(arg[10], "dn") == 0) zbcflag = 4;
  else if(strcmp(arg[10], "db") == 0) zbcflag = 5;
  else error->all(FLERR,"Illegal z-axis boundary condition command");

  rflag = 0;
//  rstep = atoi(arg[11]);
//  if (rstep > 1) rflag = 1;
}

/* ---------------------------------------------------------------------- */

FixKineticsDiffusionS::~FixKineticsDiffusionS()
{
  int i;
  for (i = 0; i < 4; i++) {
    delete [] var[i];
  }

  delete [] var;
  delete [] ivar;

  memory->destroy(diffD);
  memory->destroy(xGrid);
  memory->destroy(nuGrid);
  memory->destroy(ghost);
  memory->destroy(nuBS);
}

/* ---------------------------------------------------------------------- */

int FixKineticsDiffusionS::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsDiffusionS::init()
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 4; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix diffusion does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix diffusion is invalid style");
  }


  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  tol = input->variable->compute_equal(ivar[0]);
  q = input->variable->compute_equal(ivar[1]);
  rvol = input->variable->compute_equal(ivar[2]);
  if (rvol <= 0) lmp->error->all(FLERR,"Reactor volume cannot be equal or less than 0");
  af = input->variable->compute_equal(ivar[3]);

  //set diffusion grid size
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  nnus = bio->nnus;
  iniS = bio->iniS;
  diffCoeff = bio->diffCoeff;

  diffD = memory->create(diffD,nnus+1,"diffusion:diffD");

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi-xlo)/nx;
  stepy = (yhi-ylo)/ny;
  stepz = (zhi-zlo)/nz;
  bzhi = kinetics->bnz * stepz;

  //if (!isEuqal(stepx, stepy, stepz)) error->all(FLERR,"Grid is not cubic");

  nX = nx + 2;
  nY = ny + 2;
  nZ = nz + 2;

  nXYZ = nX*nY*nZ;

  //inlet concentration, diffusion constant
  //and maximum boundary condition conc value
  for (int i = 1; i <= nnus; i++) {
    diffD[i] = diffCoeff[i];
  }

  xGrid =  memory->create(xGrid,nXYZ,3,"diffusion:xGrid");
  nuGrid = memory->create(nuGrid,nXYZ,nnus+1,"diffusion:nuGrid");
  ghost = memory->create(ghost,nXYZ,"diffusion:ghost");
  nuBS = memory->create(nuBS,nnus+1,"diffusion:nuBS");

  //initialise grids
  double i, j, k;
  int grid = 0;
  for (k = zlo - (stepz/2); k < zhi + stepz; k += stepz) {
    for (j = ylo - (stepy/2); j < yhi + stepy; j += stepy) {
      for (i = xlo - (stepx/2); i < xhi + stepx; i += stepx) {
        xGrid[grid][0] = i;
        xGrid[grid][1] = j;
        xGrid[grid][2] = k;
        //Initialise concentration values for ghost and std grids
        for (int nu = 1; nu <= nnus; nu++) {
          if (i < xlo) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][1] * 1000;
          } else if (i > xhi) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][2] * 1000;
          } else if (j < ylo) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][3] * 1000;
          } else if (j > yhi) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][4] * 1000;
          } else if (k < zlo) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][5] * 1000;
          } else if (k > bzhi) {
            ghost[grid] = true;
            nuGrid[grid][nu] = iniS[nu][6] * 1000;
          } else {
            ghost[grid] = false;
            nuGrid[grid][nu] = iniS[nu][0] * 1000;
          }
          if (grid == 0) {
            nuBS[nu] =  iniS[nu][6] * 1000;
          }
        }
        grid++;
      }
    }
  }
}


/* ----------------------------------------------------------------------
  solve diffusion and reaction
------------------------------------------------------------------------- */

bool* FixKineticsDiffusionS::diffusion(bool *nuConv, int iter, double diffT)
{
  if (iter == 1 && kinetics->bl > 0) update_grids();

  this->diffT = diffT;
  nuS = kinetics->nuS;
  nuR = kinetics->nuR;

  //nRES = new double[nXYZ]();

  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && diffCoeff[i] != 0 && !nuConv[i]) {
      double maxS = 0;
      double *nuPrev = memory->create(nuPrev,nXYZ,"diffusion:nuPrev");
      // copy current concentrations
      for (int grid = 0; grid < nXYZ; grid++) {
        nuPrev[grid] = nuGrid[grid][i];
      }

      xbcm = iniS[i][1] * 1000;
      xbcp = iniS[i][2] * 1000;
      ybcm = iniS[i][3] * 1000;
      ybcp = iniS[i][4] * 1000;
      zbcm = iniS[i][5] * 1000;
      zbcp = iniS[i][6] * 1000;

      if(iter == 1 && strcmp(bio->nuName[i], "o2") != 0 && q >= 0 && af >= 0) compute_bulk(i);

      // solve diffusion and reaction
      for (int grid = 0; grid < nXYZ; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ix = floor(xGrid[grid][0]/stepx);
          int iy = floor(xGrid[grid][1]/stepy);
          int iz = floor(xGrid[grid][2]/stepz);

          int ind = iz * nx * ny + iy * nx + ix;

          compute_flux(diffD[i], nuGrid[grid][i], nuPrev, nuR[i][ind], grid);

          nuR[i][ind] = 0;

          if (nuGrid[grid][i] > 0) nuS[i][ind] = nuGrid[grid][i] / 1000;
          else {
            nuGrid[grid][i] = 0;
            nuS[i][ind] = 0;
          }
        } else compute_bc(nuGrid[grid][i], nuPrev, grid, nuBS[i]);

        if (maxS < nuGrid[grid][i]) maxS = nuGrid[grid][i];
      }

      bool conv = true;
      // check convergence criteria
      for (int grid = 0; grid < nXYZ; grid++) {
        if(!ghost[grid]){
          double ratio = 1000;
          if (maxS == 0) maxS = 1;

          if (rflag == 0) {
            double rate = nuGrid[grid][i]/maxS;
            double prevRate = nuPrev[grid]/maxS;
            ratio = fabs(rate - prevRate);
          }
//          else {
//            if (iter % rstep == 0) {
//              ratio = fabs(nRES[grid]);
//              nRES[grid] = 0;
//            }
//          }

          if(ratio >= tol) {
            conv = false;
            break;
          }
        }
      }
      nuConv[i] = conv;
      memory->destroy(nuPrev);
    }
  }

  //delete[] nRES;
  return nuConv;
}


void FixKineticsDiffusionS::update_grids(){
  //update grids
  bzhi = kinetics->bnz * stepz;
  nXYZ = nX * nY * (kinetics->bnz+2);

  //printf("sk = %e, grid = %i \n", sk, grid);
  for (int grid = 0; grid < nXYZ; grid++) {
    if (xGrid[grid][0] < 0 || xGrid[grid][1] < 0 || xGrid[grid][2] < 0
        || xGrid[grid][0] > xhi || xGrid[grid][1] > yhi || xGrid[grid][2] > bzhi) ghost[grid] = true;
    else ghost[grid] = false;
  }
}

/* ----------------------------------------------------------------------
  Mass balances of nutrients in the bulk liquid
------------------------------------------------------------------------- */

void FixKineticsDiffusionS::compute_bulk(int nu) {
  double sumR = 0;
  for (int i = 0; i < nx*ny*kinetics->nz; i++) sumR += nuR[nu][i];
 // printf("nuBS[%i] = %e \n", nu, nuBS[nu]);
  nuBS[nu] = nuBS[nu] + ((q/rvol) * (zbcp - nuBS[nu]) + (af/(rvol*yhi*xhi))*sumR*stepx*stepy*stepz)*update->dt * nevery;
 // printf("nuBS[%i] = %e \n", nu, nuBS[nu]);
}


/* ----------------------------------------------------------------------
  update boundary condition
------------------------------------------------------------------------- */

void FixKineticsDiffusionS::compute_bc(double &nuCell, double *nuPrev, int grid, double bulk) {
  //for nx = ny = nz = 1 grids
  //18 19 20        21 22 23       24 25 26
  //9  10 11        12 13 14       15 16 17
  //0   1  2         3  4  5        6  7  8

  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - nX;  // y direction
  int fwd = grid + nX;  // y direction
  int down = grid - nX*nY; // z direction
  int up = grid + nX*nY;  // z dirction

  // assign values to the ghost-grids according to the boundary conditions.
  // If ghostcells are Neu then take the values equal from the adjacent cells.
  // if ghostcells are dirich then take the values equal to negative of the adjacent cells.
  // if ghostcells are mixed then zlo ghost cells are nuemann, zhi ghost cells are dirichlet, other four surfaces are periodic BC.
    // Low-z surface
  if (xGrid[grid][2] < zlo && !ghost[up]) {
    //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
    if (zbcflag == 0) {
      int zhiGrid = grid + nX*nY*nz;
      nuCell = nuPrev[zhiGrid];
    } else if (zbcflag == 1) {
      nuCell = 2*zbcm - nuPrev[up];
    } else if (zbcflag == 2) {
      nuCell = nuPrev[up];
    } else if (zbcflag == 3) {
      nuCell = nuPrev[up];
    } else if (zbcflag == 4) {
      nuCell = 2*zbcm - nuPrev[up];
    }
  }
  // high-z surface
  else if (xGrid[grid][2] > zhi && !ghost[down]) {
    if (zbcflag == 0) {
      int zloGrid = grid - nX*nY*nz;
      nuCell = nuPrev[zloGrid];
    } else if (zbcflag == 1) {
      nuCell = 2*bulk - nuPrev[down];
    } else if (zbcflag == 2) {
      nuCell = 2*bulk - nuPrev[down];
    } else if (zbcflag == 3) {
      nuCell = nuPrev[down];
    } else if (zbcflag == 4) {
      nuCell = nuPrev[down];
    }
  }
  // low-y surface
  else if (xGrid[grid][1] < ylo && !ghost[fwd]) {
    if (ybcflag == 0) {
      int yhiGrid = grid + nX*ny;
      nuCell = nuPrev[yhiGrid];
    } else if (ybcflag == 1) {
      nuCell = 2*ybcm - nuPrev[fwd];
    } else if (ybcflag == 2) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == 3) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == 4) {
      nuCell = 2*ybcm - nuPrev[fwd];
    }
  }
  // high-y surface
  else if (xGrid[grid][1] > yhi && !ghost[bwd]) {
    if (ybcflag == 0) {
      int yloGrid = grid - nX*ny;
      nuCell = nuPrev[yloGrid];
    } else if (ybcflag == 1) {
      nuCell = 2*ybcp - nuPrev[bwd];
    } else if (ybcflag == 2) {
      nuCell = 2*ybcp - nuPrev[bwd];
    } else if (ybcflag == 3) {
      nuCell = nuPrev[bwd];
    } else if (ybcflag == 4) {
      nuCell = nuPrev[bwd];
    }
  }
  // low-x surface
  else if (xGrid[grid][0] < xlo && !ghost[rhs]) {
    if (xbcflag == 0) {
      int xhiGrid = grid + nx;
      nuCell = nuPrev[xhiGrid];
    } else if (xbcflag == 1) {
      nuCell = 2*xbcm - nuPrev[rhs];
    } else if (xbcflag == 2) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == 3) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == 4) {
      nuCell = 2*xbcm - nuPrev[rhs];
    }
  }
  // high-x surface
  else if (xGrid[grid][0] > xhi && !ghost[lhs]) {
    if (xbcflag == 0) {
      int xloGrid = grid - nx;
      nuCell = nuPrev[xloGrid];
    } else if (xbcflag == 1) {
      nuCell = 2*xbcp - nuPrev[lhs];
    } else if (xbcflag == 2) {
      nuCell = 2*xbcp - nuPrev[lhs];
    } else if (xbcflag == 3) {
      nuCell = nuPrev[lhs];
    } else if (xbcflag == 4) {
      nuCell = nuPrev[lhs];
    }
  }
}

/* ----------------------------------------------------------------------
  update non-ghost grids
------------------------------------------------------------------------- */

void FixKineticsDiffusionS::compute_flux(double cellDNu, double &nuCell, double *nuPrev, double rateNu, int grid) {
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - nX;  // y direction
  int fwd = grid + nX;  // y direction
  int down = grid - nX*nY; // z direction
  int up = grid + nX*nY;  // z dirction

  double jRight = cellDNu*(nuPrev[rhs] - nuPrev[grid])/stepx;
  double jLeft = cellDNu*(nuPrev[grid] - nuPrev[lhs])/stepx;
  double jX = (jRight - jLeft)/stepx;

  double jForward = cellDNu*(nuPrev[fwd] - nuPrev[grid])/stepy;
  double jBackward = cellDNu*(nuPrev[grid] - nuPrev[bwd])/stepy;
  double jY = (jForward - jBackward)/stepy;

  double jUp = cellDNu*(nuPrev[up] - nuPrev[grid])/stepz;
  double jDown = cellDNu*(nuPrev[grid] - nuPrev[down])/stepz;
  double jZ = (jUp - jDown)/stepz;

  // Adding fluxes in all the directions and the uptake rate (RHS side of the equation)
  double res = (jX + jY + jZ + rateNu) * diffT;
  //nRES[grid] += res;
  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nuCell = nuPrev[grid] + res;
}


/* ----------------------------------------------------------------------
  compare double values for equality
------------------------------------------------------------------------- */

bool FixKineticsDiffusionS::isEuqal(double a, double b, double c)
{
  double epsilon = 1e-10;
  if ((fabs(a - b) > epsilon)|| (fabs(a - b) > epsilon) || (fabs(a - b) > epsilon))
    return false;

  return true;
}

double FixKineticsDiffusionS::getMaxHeight() {
//  double minmax[6];
//  group->bounds(0,minmax);
//
//  return minmax[5];
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double **x = atom->x;
  double *r = atom->radius;
  double maxh = 0;

  for (int i=0; i<nall; i++) {
    if((x[i][2]+r[i]) > maxh) maxh = x[i][2]+r[i];
  }

  return maxh;
}


void FixKineticsDiffusionS::test(){
  //test code
  for (int i = kinetics->nz+1; i>=0; i--){
    printf(" \n ");
    for (int j = 0; j<nY; j++){
      printf("  ");
      for (int k = 0; k<nX; k++){
        int ind = i * nX * nY + j * nX + k;

//        if (!ghost[ind]) {
//          int ix = floor(xGrid[ind][0]/stepx);
//          int iy = floor(xGrid[ind][1]/stepy);
//          int iz = floor(xGrid[ind][2]/stepz);
//
//          int ind2 = iz * nx * ny + iy * nx + ix;
//
//          //printf("%i ", ind2);
//          printf("%.1e ", kinetics->nuR[1][ind2]);
//        } else {
//          printf("0 ");
//        }
        //printf("%d ", ghost[ind]);
        //printf("%i ", ind);
        //printf("%.1e ", xGrid[ind][1]);
        printf("%.1e ", nuGrid[ind][1]);
      }
    }
  }
  printf(" \n ");
}
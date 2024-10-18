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

#include "fix_density_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "grid.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "group.h"
#include "update.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDensityKokkos<DeviceType>::FixDensityKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDensity(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | MASK_MASK | BIOMASS_MASK | RMASS_MASK;
  datamask_modify = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDensityKokkos<DeviceType>::init()
{
  FixDensity::init();

  d_bitmask = typename AT::t_int_1d("nufeb/density/kk:bitmask", group->ngroup);
  // copy group bitmask to device
  Kokkos::View<int *, Kokkos::HostSpace,
	       Kokkos::MemoryTraits<Kokkos::Unmanaged> > h_bitmask(group->bitmask, group->ngroup);
  Kokkos::deep_copy(d_bitmask, h_bitmask);
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void FixDensityKokkos<DeviceType>::post_physics_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}
/* ---------------------------------------------------------------------- */

struct AtomicD {
  enum {value = Kokkos::Atomic|Kokkos::Unmanaged};
};

template<class DeviceType>
void FixDensityKokkos<DeviceType>::compute()
{
  atomKK->sync(execution_space,datamask_read);

  for (int i = 0; i < 3; i++) {
    boxlo[i] = domain->boxlo[i];
    sublo[i] = domain->sublo[i];
    subhi[i] = domain->subhi[i];
    grid_sublo[i] = grid->sublo[i];
    grid_subbox[i] = grid->subbox[i];
  }
  ngroup = group->ngroup;
  cell_size = grid->cell_size;
  vol = cell_size * cell_size * cell_size;

  d_mask = atomKK->k_mask.view<DeviceType>();
  d_x = atomKK->k_x.view<DeviceType>();
  d_biomass = atomKK->k_biomass.template view<DeviceType>();
  d_rmass = atomKK->k_rmass.template view<DeviceType>();
  d_dens = gridKK->k_dens.template view<DeviceType>();

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDensityInitialTag>(0, grid->ncells), f);

  DeviceType().fence();

  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDensityComputeTag>(0, atomKK->nlocal+atomKK->nghost), f);
  DeviceType().fence();
  copymode = 0;

  gridKK->modified(execution_space, DENS_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDensityKokkos<DeviceType>::Functor::Functor(FixDensityKokkos<DeviceType> *ptr):
  ngroup(ptr->ngroup), d_dens(ptr->d_dens),
  cell_size(ptr->cell_size), vol(ptr->vol),
  d_mask(ptr->d_mask), d_x(ptr->d_x),
  d_biomass(ptr->d_biomass), d_rmass(ptr->d_rmass),
  d_bitmask(ptr->d_bitmask)
{
  for (int i = 0; i < 3; i++) {
    boxlo[i] = ptr->boxlo[i];
    sublo[i] = ptr->sublo[i];
    subhi[i] = ptr->subhi[i];
    grid_sublo[i] = ptr->grid_sublo[i];
    grid_subbox[i] = ptr->grid_subbox[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDensityKokkos<DeviceType>::Functor::operator()(FixDensityInitialTag, int i) const
{
  for (int igroup = 0; igroup < ngroup; igroup++) {
    d_dens(igroup,i) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDensityKokkos<DeviceType>::Functor::operator()(FixDensityComputeTag, int i) const
{
  Kokkos::View<double**,
	      typename DAT::t_float_2d::array_layout,DeviceType,
	      Kokkos::MemoryTraits<AtomicD::value> > a_dens = d_dens;

  // including ghost atoms because there can be atoms that moved inside the
  //   sub-domain and were not yet exchanged
  // forward communication garantees that we have the latest ghost positions
  //   which were updated during initial integrate
  if (d_x(i,0) >= sublo[0] && d_x(i,0) < subhi[0] &&
     d_x(i,1) >= sublo[1] && d_x(i,1) < subhi[1] &&
     d_x(i,2) >= sublo[2] && d_x(i,2) < subhi[2]) {
   // can't do:
   // int cell = grid->cell(d_x[i]);
   // hence copying the code here
   int c[3];
   const double small = 1e-12;
   for (int j = 0; j < 3; j++) {
     c[j] = static_cast<int>((d_x(i,j) - boxlo[j]) /
			     cell_size + small) - grid_sublo[j];
   }
   int cell = c[0] + c[1] * grid_subbox[0] +
     c[2] * grid_subbox[0] * grid_subbox[1];

   double d = d_biomass(i) * d_rmass(i) / vol;

   a_dens(0,cell) += d;
   for (int igroup = 0; igroup < ngroup; igroup++)
     if (d_mask(i) & d_bitmask(igroup))
       a_dens(igroup,cell) += d;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDensityKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDensityKokkos<LMPHostType>;
#endif
}

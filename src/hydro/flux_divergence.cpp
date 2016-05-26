//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file flux_divergence.cpp
//  \brief Applies divergence of the fluxes, including geometric "source terms" added
//         by a function implemented in each Coordinate class.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../mesh.hpp"
#include "srcterms/srcterms.hpp"
#include "../bvals/bvals.hpp"
#include "../reconstruct/reconstruction.hpp"

// this class header
#include "hydro.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//--------------------------------------------------------------------------------------
//! \fn  void Hydro::FluxDivergence
//  \brief Integrate the conservative variables using the calculated fluxes

void Hydro::FluxDivergence(MeshBlock *pmb,AthenaArray<Real> &u,
  AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &bcc, const int step)
{
  AthenaArray<Real> &x1flux=flux[x1face];
  AthenaArray<Real> &x2flux=flux[x2face];
  AthenaArray<Real> &x3flux=flux[x3face];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dt;
  if (step == 1) {
    dt = 0.5*(pmb->pmy_mesh->dt);
  } else {
    dt = (pmb->pmy_mesh->dt);
  }

  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
  x1area.InitWithShallowSlice(x1face_area_,2,tid,1);
  x2area.InitWithShallowSlice(x2face_area_,2,tid,1);
  x2area_p1.InitWithShallowSlice(x2face_area_p1_,2,tid,1);
  x3area.InitWithShallowSlice(x3face_area_,2,tid,1);
  x3area_p1.InitWithShallowSlice(x3face_area_p1_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);
  dflx.InitWithShallowSlice(flx_,2,tid,1);

#pragma omp for schedule(static)
  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for (int n=0; n<NHYDRO; ++n) {
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int n=0; n<NHYDRO; ++n) {
        for (int i=is; i<=ie; ++i) {
          u(n,k,j,i) -= dt*dflx(n,i)/vol(i);
        }
      }
    }
  }

} // end of omp parallel region


  // add coordinate (geometric) source terms
  pmb->pcoord->CoordSrcTerms(dt,pmb->phydro->flux,w,bcc,u);

  return;
}

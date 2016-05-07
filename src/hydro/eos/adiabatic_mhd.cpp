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
//! \file adiabatic_mhd.cpp
//  \brief implements functions in class HydroEqnOfState for adiabatic MHD
//======================================================================================

// C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN

// Athena++ headers
#include "../hydro.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../field/field.hpp"
#include "../../coordinates/coordinates.hpp"

// this class header
#include "eos.hpp"

// HydroEqnOfState constructor

HydroEqnOfState::HydroEqnOfState(Hydro *pf, ParameterInput *pin)
{
  pmy_hydro_ = pf;
  gamma_ = pin->GetReal("hydro","gamma");
  density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
}

// destructor

HydroEqnOfState::~HydroEqnOfState()
{
}

//--------------------------------------------------------------------------------------
// \!fn void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//    const AthenaArray<Real> &prim_old, const FaceField &b,
//    AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//    int is, int ie, int js, int je, int ks, int ke);
// \brief For the Hydro, converts conserved into primitive variables in adiabatic MHD.
//  For the Field, computes cell-centered from face-centered magnetic field.

void HydroEqnOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
	    const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
    AthenaArray<Real> &bcc, Coordinates *pco,
    int is, int ie, int js, int je, int ks, int ke)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Real gm1 = GetGamma() - 1.0;

  pmb->pfield->CalculateCellCenteredField(b,bcc,pco,is,ie,js,je,ks,ke);

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
//#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IVX,k,j,i);
      Real& u_m2 = cons(IVY,k,j,i);
      Real& u_m3 = cons(IVZ,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      Real& w_d  = prim(IDN,k,j,i);
      Real& w_vx = prim(IVX,k,j,i);
      Real& w_vy = prim(IVY,k,j,i);
      Real& w_vz = prim(IVZ,k,j,i);
      Real& w_p  = prim(IEN,k,j,i);

      // apply density floor, without changing momentum or energy
      u_d = (u_d > density_floor_) ?  u_d : density_floor_;
      w_d = u_d;

      Real di = 1.0/u_d;
      w_vx = u_m1*di;
      w_vy = u_m2*di;
      w_vz = u_m3*di;

      const Real& bcc1 = bcc(IB1,k,j,i);
      const Real& bcc2 = bcc(IB2,k,j,i);
      const Real& bcc3 = bcc(IB3,k,j,i);

      Real pb = 0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
      Real ke = 0.5*di*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
      w_p = gm1*(u_e - ke - pb);

      // apply pressure floor, correct total energy
      u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1) + ke + pb);
      w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
    }
  }}
}

  return;
}

//--------------------------------------------------------------------------------------
// \!fn void HydroEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables
//        Note that this function assumes cell-centered fields are already calculated

void HydroEqnOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Real igm1 = 1.0/(GetGamma() - 1.0);

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
//#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& u_d  = cons(IDN,k,j,i);
      Real& u_m1 = cons(IM1,k,j,i);
      Real& u_m2 = cons(IM2,k,j,i);
      Real& u_m3 = cons(IM3,k,j,i);
      Real& u_e  = cons(IEN,k,j,i);

      const Real& w_d  = prim(IDN,k,j,i);
      const Real& w_vx = prim(IVX,k,j,i);
      const Real& w_vy = prim(IVY,k,j,i);
      const Real& w_vz = prim(IVZ,k,j,i);
      const Real& w_p  = prim(IEN,k,j,i);

      const Real& bcc1 = bc(IB1,k,j,i);
      const Real& bcc2 = bc(IB2,k,j,i);
      const Real& bcc3 = bc(IB3,k,j,i);

      u_d = w_d;
      u_m1 = w_vx*w_d;
      u_m2 = w_vy*w_d;
      u_m3 = w_vz*w_d;
      u_e = w_p*igm1 + 0.5*(w_d*(SQR(w_vx) + SQR(w_vy) + SQR(w_vz))
            + (SQR(bcc1) + SQR(bcc2) + SQR(bcc3)));
    }
  }}
}
  return;
}

//--------------------------------------------------------------------------------------
// \!fn Real HydroEqnOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables

Real HydroEqnOfState::SoundSpeed(const Real prim[NHYDRO])
{
  return sqrt(GetGamma()*prim[IEN]/prim[IDN]);
}

//--------------------------------------------------------------------------------------
// \!fn Real HydroEqnOfState::FastMagnetosonicSpeed(const Real prim[], const Real bx)
// \brief returns fast magnetosonic speed given vector of primitive variables
// Note the formula for (C_f)^2 is positive definite, so this func never returns a NaN 

Real HydroEqnOfState::FastMagnetosonicSpeed(const Real prim[(NWAVE)], const Real bx)
{
  Real asq = GetGamma()*prim[IEN]/prim[IDN];
  Real vaxsq = bx*bx/prim[IDN];
  Real ct2 = (prim[IBY]*prim[IBY] + prim[IBZ]*prim[IBZ])/prim[IDN];
  Real qsq = vaxsq + ct2 + asq;
  Real tmp = vaxsq + ct2 - asq;
  return sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ct2)));
}

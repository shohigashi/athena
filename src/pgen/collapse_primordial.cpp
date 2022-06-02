//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file collapse_primordial.cpp
//! \brief Problem generator for collapse of a Bonnor-Ebert like sphere with AMR or SMR

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()
#include "../fft/athena_fft.hpp"

#if SELF_GRAVITY_ENABLED != 2
#error "This problem generator requires Multigrid gravity solver."
#endif


namespace {


// dimension-less constants
constexpr Real four_pi_G = 1.0;
constexpr Real rc = 6.45; // the BE radius in the normalized unit system
constexpr Real rcsq = 26.0 / 3.0;      // the parameter of the BE profile
constexpr Real bemass = 197.561;       // the total mass of the critical BE sphere

// dimensional constants
constexpr Real pi   = M_PI;
constexpr Real cs10 = 1.9e4;        // sound speed at 10K, cm / s 
constexpr Real cs200 = 1.214e5;        // sound speed at 200K, cm / s 
constexpr Real msun = 1.9891e33;    // solar mass, g
constexpr Real msune3 = 1.9891e36;    // 1e3 solar mass, g
constexpr Real pc   = 3.0857000e18; // parsec, cm
constexpr Real au   = 1.4959787e13; // astronomical unit, cm
constexpr Real yr   = 3.15569e7;    // year, s
constexpr Real G    = 6.67259e-8;   // gravitational constant, dyn cm^2 g^-2

// units in cgs
Real m0, v0, t0, l0, rho0;

// parameters and derivatives
//Real mass, temp, f, rhocrit;
Real mass, temp, f, rhocrit, b0, geff; // b-field is added
// Real angle; // reserved

// AMR parameter
Real njeans; // Real is used intentionally 

}

// Mask the density outside the initial sphere
void SourceMask(AthenaArray<Real> &src, int is, int ie, int js, int je,
                int ks, int ke, const MGCoordinates &coord) {
  const Real rc2 = rc*rc;
  for (int k=ks; k<=ke; ++k) {
    Real z = coord.x3v(k);
    for (int j=js; j<=je; ++j) {
      Real y = coord.x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x = coord.x1v(i);
        Real r2 = x*x + y*y + z*z;
        if (r2 > rc2)
          src(k, j, i) = 0.0;
      }
    }
  }
  return;
}


// AMR refinement condition
int JeansCondition(MeshBlock *pmb) {
  Real njmin = 1e300;
  const Real dx = pmb->pcoord->dx1f(0); // assuming uniform cubic cells
  /* Jeans criterion is the same not dependent on the existence of B-fields for comparison with the model without B-fields.
   // Sho Higashi 02/06/2022
  if (MAGNETIC_FIELDS_ENABLED) {
    if (NON_BAROTROPIC_EOS) {
      const Real gamma = pmb->peos->GetGamma();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real v = std::sqrt(gamma * pmb->phydro->w(IPR,k,j,i)
                             +(SQR(pmb->pfield->bcc(IB1,k,j,i))
                             + SQR(pmb->pfield->bcc(IB2,k,j,i))
                             + SQR(pmb->pfield->bcc(IB3,k,j,i)))
                             / pmb->phydro->w(IDN,k,j,i));
            Real nj = fac * v / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
    } else {
      const Real cs = pmb->peos->GetIsoSoundSpeed();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real v = cs + std::sqrt((SQR(pmb->pfield->bcc(IB1,k,j,i))
                                   + SQR(pmb->pfield->bcc(IB2,k,j,i))
                                   + SQR(pmb->pfield->bcc(IB3,k,j,i)))
                                   / pmb->phydro->w(IDN,k,j,i));
            Real nj = fac * v / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
    }
  } else {*/
    if (NON_BAROTROPIC_EOS) {
      const Real gamma = pmb->peos->GetGamma();
      const Real fac = 2.0 * pi * std::sqrt(gamma) / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real nj = fac * std::sqrt(pmb->phydro->w(IPR,k,j,i))
                              / pmb->phydro->w(IDN,k,j,i);
            njmin = std::min(njmin, nj);
          }
        }
      }
    } else {
      const Real cs = pmb->peos->GetIsoSoundSpeed();
      const Real fac = 2.0 * pi / dx;
      for (int k = pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
        for (int j = pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
          for (int i = pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
            Real nj = fac * cs / std::sqrt(pmb->phydro->w(IDN,k,j,i));
            njmin = std::min(njmin, nj);
          }
        }
      }
    }
 // }
  if (njmin < njeans)
    return 1;
  if (njmin > njeans * 2.5)
    return -1;
  return 0;
}

// Approximated Bonnor-Ebert profile
// Tomida 2011, PhD Thesis
Real BEProfile(Real r) {
  return std::pow(1.0+r*r/rcsq, -1.5);
}


void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  const Real gm1 = pmb->peos->GetGamma() - 1.0;
  const Real igm1 = 1.0 / gm1;
  const Real geff1 = geff - 1.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real ke = 0.5 / cons(IDN,k,j,i)
                  * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
          Real me = 0.5*(SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)));
//          Real te = igm1 * cons(IDN,k,j,i)
//                  * std::max(1.0, std::pow(cons(IDN,k,j,i)/rhocrit, gm1));
          Real te = igm1 * cons(IDN,k,j,i)
                  * std::max(std::pow(cons(IDN,k,j,i)/rhocrit,geff1), std::pow(cons(IDN,k,j,i)/rhocrit, gm1))
                  / std::pow(1.0/rhocrit,geff1);
          cons(IEN,k,j,i) = te + ke + me;
        }
      }
    }
  } else {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real ke = 0.5 / cons(IDN,k,j,i)
                  * (SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)));
//          Real te = igm1 * cons(IDN,k,j,i)
//                  * std::max(1.0, std::pow(cons(IDN,k,j,i)/rhocrit, gm1));
          Real te = igm1 * cons(IDN,k,j,i)
                  * std::max(std::pow(cons(IDN,k,j,i)/rhocrit,geff1), std::pow(cons(IDN,k,j,i)/rhocrit, gm1))
                  / std::pow(1.0/rhocrit,geff1);
          cons(IEN,k,j,i) = te + ke;
        }
      }
    }
  }

  return;
}




//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  SetFourPiG(four_pi_G); // 4piG = 1.0
  mass = pin->GetReal("problem", "mass");
  temp = pin->GetReal("problem", "temperature");
  f = pin->GetReal("problem", "f"); // Density enhancement factor; f = 1 is critical
  geff = pin->GetReal("problem","geff");
  m0 = mass * msune3 / (bemass*f);
  //v0 = cs10 * std::sqrt(temp/10.0);
  v0 = cs200 * std::sqrt(temp/200.0);
  rho0 = (v0*v0*v0*v0*v0*v0) / (m0*m0) /(64.0*pi*pi*pi*G*G*G);
  t0 = 1.0/std::sqrt(4.0*pi*G*rho0);
  l0 = v0 * t0;
  rhocrit = pin->GetReal("problem", "rhocrit") / rho0;
  Real tff = sqrt(3.0/8.0)*pi;
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  if (MAGNETIC_FIELDS_ENABLED){
    b0  = pin->GetReal("problem", "b0");
  }else{
    b0 = 0.0;
  }

  if (Globals::my_rank == 0 && ncycle == 0) {
    std::cout << std::endl
      << "---  Dimensional parameters of the simulation  ---" << std::endl
      << "Total mass          : " << mass*1e3      << " \t\t[Msun]" << std::endl
      << "Initial temperature : " << temp      << " \t\t[K]" << std::endl
      << "Sound speed         : " << v0        << " \t\t[cm s^-1]" << std::endl
      << "Central density     : " << rho0      << " \t[g cm^-3]" << std::endl
      << "Cloud radius        : " << rc*l0/pc  << " \t\t[pc]" << std::endl
      << "Free fall time      : " << tff*t0/yr << " \t\t[yr]" << std::endl
      << "Density Enhancement : " << f         << std::endl
      << "Magnetic field      : " << b0*std::sqrt(4*M_PI)*std::sqrt(rho0*SQR(v0)) << " \t\t[G]" << std::endl
      << std::endl
      << "---   Normalization Units of the simulation    ---" << std::endl
      << "Mass                : " << m0        << " \t[g]" << std::endl
      << "Mass                : " << m0/msune3   << " \t[Msun]" << std::endl
      << "Length              : " << l0        << " \t[cm]" << std::endl
      << "Length              : " << l0/au     << " \t\t[au]" << std::endl
      << "Length              : " << l0/pc     << " \t[pc]" << std::endl
      << "Time                : " << t0        << " \t[s]" << std::endl
      << "Time                : " << t0/yr     << " \t\t[yr]" << std::endl
      << "Velocity            : " << v0        << " \t\t[cm s^-1]" << std::endl
      << "Density             : " << rho0      << " \t[g cm^-3]" << std::endl
      << "Magnetic field      : " << b0*std::sqrt(4*M_PI)*std::sqrt(rho0*SQR(v0)) << " \t\t[G]" << std::endl
      << std::endl
      << "--- Dimensionless parameters of the simulation ---" << std::endl
      << "Total mass          : " << bemass*f  << std::endl
      << "Sound speed at " << temp << " K : "  << 1.0 << std::endl
      << "Central density     : " << 1.0       << std::endl
      << "Cloud radius        : " << rc        << std::endl
      << "Free fall time      : " << tff       << std::endl
      << std::endl;
  }

  EnrollUserMGGravitySourceMaskFunction(SourceMask);

  if (NON_BAROTROPIC_EOS)
    EnrollUserExplicitSourceFunction(Cooling);

  if (adaptive) {
    njeans = pin->GetReal("problem","njeans");
    EnrollUserRefinementCondition(JeansCondition);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real igm1 = 1.0 / (peos->GetGamma() - 1.0);
  for (int k=ks; k<=ke; ++k) {
    Real z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real y = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x = pcoord->x1v(i);
        Real r = std::sqrt(SQR(x) + SQR(y) + SQR(z));
        r = std::min(r, rc); // pressure confinement - constant beyond the cloud radius
        phydro->u(IDN,k,j,i) = BEProfile(r);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS)
          phydro->u(IEN,k,j,i) = igm1 * phydro->u(IDN,k,j,i); // c_s = 1
      }
    }
  }
  
  if (MAGNETIC_FIELDS_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
        pfield->b.x1f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
        pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
        pfield->b.x3f(k,j,i) = b0;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*SQR(b0);
        }
      }
    }
  }
  // Determine locations of initial values
  //int USE_TURB_FROM_HDF5 = pin->GetOrAddBoolean("problem","use_turb_from_hdf5",false);
  //if (USE_TURB_FROM_HDF5){
  
  /*
  if (USE_TURB_FROM_HDF5){
      std::string input_filename = pin->GetString("problem", "input_filename");
      std::string dataset_v1 = pin->GetString("problem", "dataset_v1");
      std::string dataset_v2 = pin->GetString("problem", "dataset_v2");
      std::string dataset_v3 = pin->GetString("problem", "dataset_v3");
      int turb_dim = pin->GetInteger("problem","turb_dim");
 */   /*
      int start_field_file[4] = {gid,0,0,0};
      //int count_field_file[4] = {1,turb_dim,turb_dim,turb_dim};
      int count_field_file[4] = {1,block_size.nx3,block_size.nx2,block_size.nx1};
      //int start_field_mem[3] = {0,0,0};
      int start_field_mem[3] = {ks,js,is};
      //int count_field_mem[3] = {turb_dim,turb_dim,turb_dim};
      int count_field_mem[3] = {block_size.nx3,block_size.nx2,block_size.nx1};
      if (Globals::my_rank == 0) {
        for (int n=0;n<3;++n){
        std::cout << "block_size:" << count_field_file[n] << std::endl; 
        std::cout << "ks,js,is:" <<  start_field_mem[n] << std::endl; 
        }
      }
      // Set field array selections
      AthenaArray<Real> turb_x;
      //turb_x.NewAthenaArray(ke-ks,je-js,ie-is);
      turb_x.NewAthenaArray(turb_dim,turb_dim,turb_dim);
      HDF5ReadRealArray(input_filename.c_str(), dataset_v1.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, turb_x,true);//read x-velocity
      AthenaArray<Real> turb_y;
      //turb_y.NewAthenaArray(ke-ks,je-js,ie-is);
      turb_y.NewAthenaArray(turb_dim,turb_dim,turb_dim);
      HDF5ReadRealArray(input_filename.c_str(), dataset_v2.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, turb_y,true);//read y-velocity
      AthenaArray<Real> turb_z;
      //turb_z.NewAthenaArray(ke-ks,je-js,ie-is);
      turb_z.NewAthenaArray(turb_dim,turb_dim,turb_dim);
      HDF5ReadRealArray(input_filename.c_str(), dataset_v3.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, turb_z,true);//read y-velocity
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
              phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*turb_x(k,j,i);
              phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*turb_y(k,j,i);
              phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*turb_z(k,j,i);
          }
        }
      }
        turb_x.DeleteAthenaArray();
        turb_y.DeleteAthenaArray();
        turb_z.DeleteAthenaArray();
    }else{
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IM1,k,j,i) = 0.0;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = 0.0;
          }
        }
      }
    }
        */
}
//========================================================================================
/* User defined output */
//========================================================================================

/*
//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief
//========================================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
    AllocateUserOutputVariables(2); // temperature, plasma beta
    return;
}
//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief
//========================================================================================
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
    for(int k=ks; k<=ke; k++) {
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie; i++) {
           Real pmag = 0.5*(SQR(pfield->bcc(IB1,k,j,i))
                     +SQR(pfield->bcc(IB2,k,j,i))
                     +SQR(pfield->bcc(IB3,k,j,i)));
           user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);
           user_out_var(1,k,j,i) = phydro->w(IPR,k,j,i)/pmag;
        }
      }
    }
}

*/

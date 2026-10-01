//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_blockdt.cpp
//! \brief computes timestep using CFL condition on a MEshBlock

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // abs(), sqrt()
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/implicit/radiation_implicit.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../scalars/scalars.hpp"
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn void Hydro::NewBlockTimeStep()
//! \brief calculate the minimum timestep within a MeshBlock

void Hydro::NewBlockTimeStep() {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &w = pmb->phydro->w;
  // hyperbolic timestep constraint in each (x1-slice) cell along coordinate direction:
  AthenaArray<Real> &dt1 = dt1_, &dt2 = dt2_, &dt3 = dt3_;  // (x1 slices)
  Real wi[NWAVE];

  Real real_max = std::numeric_limits<Real>::max();
  Real min_dt = real_max;
  // Note, "dt_hyperbolic" currently refers to the dt limit imposed by evoluiton of the
  // ideal hydro or MHD fluid by the main integrator (even if not strictly hyperbolic)
  Real min_dt_hyperbolic  = real_max;
  // TODO(felker): consider renaming dt_hyperbolic after general execution model is
  // implemented and flexibility from #247 (zero fluid configurations) is
  // addressed. dt_hydro, dt_main (inaccurate since "dt" is actually main), dt_MHD?
  Real min_dt_parabolic  = real_max;
  Real min_dt_user  = real_max;

  // Diagnostic target: first active polar cell in gid 0.
  const bool debug_block = (pmb->gid == 0);
  polar_debug_.cfl_valid = false;

  Real cspeed = 0.0;
  if(NR_RADIATION_ENABLED)
    cspeed = pmb->pnrrad->reduced_c;
  if(CR_ENABLED)
    cspeed = std::max(cspeed,pmb->pcr->vmax);

  // TODO(felker): skip this next loop if pm->fluid_setup == FluidFormulation::disabled
  FluidFormulation fluid_status = pmb->pmy_mesh->fluid_setup;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // Obtain cell widths
      pmb->pcoord->CenterWidth1(k, j, is, ie, dt1);
      pmb->pcoord->CenterWidth2(k, j, is, ie, dt2);
      pmb->pcoord->CenterWidth3(k, j, is, ie, dt3);

      if (debug_block && k == ks && j == js) {
        polar_debug_.center_width1 = dt1(is);
        polar_debug_.center_width2 = dt2(is);
      }

      // Newtonian case: divide cell widths by maximum characteristic speed
      if (!RELATIVISTIC_DYNAMICS) {
#pragma ivdep
        for (int i=is; i<=ie; ++i) {
          wi[IDN] = w(IDN,k,j,i);
          wi[IVX] = w(IVX,k,j,i);
          wi[IVY] = w(IVY,k,j,i);
          wi[IVZ] = w(IVZ,k,j,i);
          if (NON_BAROTROPIC_EOS) wi[IPR] = w(IPR,k,j,i);
          if (fluid_status == FluidFormulation::evolve) {
            if (MAGNETIC_FIELDS_ENABLED) {
              AthenaArray<Real> &bcc = pmb->pfield->bcc, &b_x1f = pmb->pfield->b.x1f,
                              &b_x2f = pmb->pfield->b.x2f, &b_x3f = pmb->pfield->b.x3f;
              Real bx = bcc(IB1,k,j,i) + std::abs(b_x1f(k,j,i) - bcc(IB1,k,j,i));
              wi[IBY] = bcc(IB2,k,j,i);
              wi[IBZ] = bcc(IB3,k,j,i);
              Real cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
              Real speed = std::max(cspeed,(std::abs(wi[IVX]) + cf));
              dt1(i) /= (speed);

              wi[IBY] = bcc(IB3,k,j,i);
              wi[IBZ] = bcc(IB1,k,j,i);
              bx = bcc(IB2,k,j,i) + std::abs(b_x2f(k,j,i) - bcc(IB2,k,j,i));
              cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
              speed = std::max(cspeed,(std::abs(wi[IVY]) + cf));
              dt2(i) /= (speed);

              wi[IBY] = bcc(IB1,k,j,i);
              wi[IBZ] = bcc(IB2,k,j,i);
              bx = bcc(IB3,k,j,i) + std::abs(b_x3f(k,j,i) - bcc(IB3,k,j,i));
              cf = pmb->peos->FastMagnetosonicSpeed(wi,bx);
              speed = std::max(cspeed,(std::abs(wi[IVZ]) + cf));
              dt3(i) /= (speed);
            } else {
#if EOS_SCALAR_INPUT_ENABLED
              Real r_cell[(NSCALARS > 0) ? NSCALARS : 1];
              for (int n=0; n<NSCALARS; ++n) {
                r_cell[n] = pmb->pscalars->r(n,k,j,i);
              }
              Real cs = std::sqrt(pmb->peos->AsqFromRhoP(wi[IDN], wi[IPR], r_cell));
#else
              Real cs = pmb->peos->SoundSpeed(wi);
#endif
              Real speed1 = std::max(cspeed, (std::abs(wi[IVX]) + cs));
              Real speed2 = std::max(cspeed, (std::abs(wi[IVY]) + cs));
              Real speed3 = std::max(cspeed, (std::abs(wi[IVZ]) + cs));
              dt1(i) /= (speed1);
              dt2(i) /= (speed2);
              dt3(i) /= (speed3);
            }
          } else { // FluidFormulation::background or disabled. Assume scalar advection:
            dt1(i) /= (std::abs(wi[IVX]));
            dt2(i) /= (std::abs(wi[IVY]));
            dt3(i) /= (std::abs(wi[IVZ]));
          }
        }
      }

      // SR case: do nothing (assume maximum characteristic is c = 1)
      // GR case: divide cell widths by coordinate speed of light (not necessarily unity)
      if (GENERAL_RELATIVITY) {
        pmb->pcoord->CellMetric(k, j, is, ie, g_, gi_);
        for (int i=is; i<=ie; ++i) {
          Real speed1 = -(std::sqrt(SQR(gi_(I01,i)) - gi_(I00,i) * gi_(I11,i))
              + std::abs(gi_(I01,i))) / gi_(I00,i);
          Real speed2 = -(std::sqrt(SQR(gi_(I02,i)) - gi_(I00,i) * gi_(I22,i))
              + std::abs(gi_(I02,i))) / gi_(I00,i);
          Real speed3 = -(std::sqrt(SQR(gi_(I03,i)) - gi_(I00,i) * gi_(I33,i))
              + std::abs(gi_(I03,i))) / gi_(I00,i);
          dt1(i) /= speed1;
          dt2(i) /= speed2;
          dt3(i) /= speed3;

          if (debug_block && k == ks && j == js && i == is) {
            Coordinates *pco = pmb->pcoord;
            const Real vol = pco->GetCellVolume(k, j, i);
            const Real a1_l = pco->GetFace1Area(k, j, i);
            const Real a1_r = pco->GetFace1Area(k, j, i+1);
            const Real a2_l = pco->GetFace2Area(k, j, i);
            const Real a2_r = pco->GetFace2Area(k, j+1, i);
            const Real rate1 = std::max(a1_l, a1_r)*speed1/vol;
            const Real rate2 = std::max(a2_l, a2_r)*speed2/vol;

            polar_debug_.cfl_valid = true;
            polar_debug_.cfl_cycle = pmb->pmy_mesh->ncycle;
            polar_debug_.cfl_time = pmb->pmy_mesh->time;
            polar_debug_.r = pco->x1v(i);
            polar_debug_.theta = pco->x2v(j);
            polar_debug_.dtheta = pco->dx2f(j);
            polar_debug_.g00 = g_(I00,i);
            polar_debug_.g02 = g_(I02,i);
            polar_debug_.g22 = g_(I22,i);
            polar_debug_.gi00 = gi_(I00,i);
            polar_debug_.gi02 = gi_(I02,i);
            polar_debug_.gi22 = gi_(I22,i);
            polar_debug_.speed1 = speed1;
            polar_debug_.speed2 = speed2;
            polar_debug_.dt1_raw = dt1(i);
            polar_debug_.dt2_raw = dt2(i);
            polar_debug_.rho = w(IDN,k,j,i);
            polar_debug_.press = w(IPR,k,j,i);
            polar_debug_.u1 = w(IVX,k,j,i);
            polar_debug_.u2 = w(IVY,k,j,i);
            polar_debug_.u3 = w(IVZ,k,j,i);
            polar_debug_.volume = vol;
            polar_debug_.area1_l = a1_l;
            polar_debug_.area1_r = a1_r;
            polar_debug_.area2_l = a2_l;
            polar_debug_.area2_r = a2_r;
            polar_debug_.fv_dt1_raw = 1.0/rate1;
            polar_debug_.fv_dt2_raw = 1.0/rate2;
            polar_debug_.fv_dt_multid_scaled =
                pmb->pmy_mesh->cfl_number/(rate1 + rate2);
          }
        }
      }

      // compute minimum of (v1 +/- C)
      for (int i=is; i<=ie; ++i) {
        const Real& dt_1 = dt1(i);
        min_dt_hyperbolic = std::min(min_dt_hyperbolic, dt_1);
      }

      // if grid is 2D/3D, compute minimum of (v2 +/- C)
      if (pmb->block_size.nx2 > 1) {
        for (int i=is; i<=ie; ++i) {
          const Real& dt_2 = dt2(i);
          min_dt_hyperbolic = std::min(min_dt_hyperbolic, dt_2);
        }
      }

      // if grid is 3D, compute minimum of (v3 +/- C)
      if (pmb->block_size.nx3 > 1) {
        for (int i=is; i<=ie; ++i) {
          const Real& dt_3 = dt3(i);
          min_dt_hyperbolic = std::min(min_dt_hyperbolic, dt_3);
        }
      }
    }
  }

  if (polar_debug_.cfl_valid) {
    polar_debug_.block_dt_raw = min_dt_hyperbolic;
  }

  // calculate the timestep limited by the diffusion processes
  if (hdif.hydro_diffusion_defined) {
    Real min_dt_vis, min_dt_cnd;
    hdif.NewDiffusionDt(min_dt_vis, min_dt_cnd);
    min_dt_parabolic = std::min(min_dt_parabolic, min_dt_vis);
    min_dt_parabolic = std::min(min_dt_parabolic, min_dt_cnd);
  } // hydro diffusion

  if (MAGNETIC_FIELDS_ENABLED &&
      pmb->pfield->fdif.field_diffusion_defined) {
    Real min_dt_oa, min_dt_hall;
    pmb->pfield->fdif.NewDiffusionDt(min_dt_oa, min_dt_hall);
    min_dt_parabolic = std::min(min_dt_parabolic, min_dt_oa);
    // Hall effect is dispersive, not diffusive:
    min_dt_hyperbolic = std::min(min_dt_hyperbolic, min_dt_hall);
  } // field diffusion

  if (NSCALARS > 0 && pmb->pscalars->scalar_diffusion_defined) {
    Real min_dt_scalar_diff = pmb->pscalars->NewDiffusionDt();
    min_dt_parabolic = std::min(min_dt_parabolic, min_dt_scalar_diff);
  } // passive scalar diffusion

  min_dt_hyperbolic *= pmb->pmy_mesh->cfl_number;
  if (polar_debug_.cfl_valid) {
    polar_debug_.block_dt_scaled = min_dt_hyperbolic;
  }
  // scale the theoretical stability limit by a safety factor = the hyperbolic CFL limit
  // (user-selected or automaticlaly enforced). May add independent parameter "cfl_diff"
  // in the future (with default = cfl_number).
  min_dt_parabolic *= pmb->pmy_mesh->cfl_number;

  if (IM_RADIATION_ENABLED) {
    min_dt_hyperbolic *= pmb->pmy_mesh->pimrad->cfl_rad;
    min_dt_parabolic  *= pmb->pmy_mesh->pimrad->cfl_rad;
  }

  // For orbital advection, give a restriction on dt_hyperbolic.
  if (pmb->porb->orbital_advection_active) {
    Real min_dt_orb = pmb->porb->NewOrbitalAdvectionDt();
    min_dt_hyperbolic = std::min(min_dt_hyperbolic, min_dt_orb);
  }

  // set main integrator timestep as the minimum of the appropriate timestep constraints:
  // hyperbolic: (skip if fluid is nonexistent or frozen)
  min_dt = std::min(min_dt, min_dt_hyperbolic);
  // user:
  if (UserTimeStep_ != nullptr) {
    min_dt_user = UserTimeStep_(pmb);
    min_dt = std::min(min_dt, min_dt_user);
  }
  // parabolic:
  // STS handles parabolic terms -> then take the smaller of hyperbolic or user timestep
  if (!STS_ENABLED) {
    // otherwise, take the smallest of the hyperbolic, parabolic, user timesteps
    min_dt = std::min(min_dt, min_dt_parabolic);
  }
  pmb->new_block_dt_ = min_dt;
  pmb->new_block_dt_hyperbolic_ = min_dt_hyperbolic;
  pmb->new_block_dt_parabolic_ = min_dt_parabolic;
  pmb->new_block_dt_user_ = min_dt_user;

  return;
}

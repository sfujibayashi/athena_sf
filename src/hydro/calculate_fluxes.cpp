//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_fluxes.cpp
//! \brief Calculate hydro/MHD fluxes

// C headers
#include <iostream>   // endl
#include <sstream>    // stringstream

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes
//! \brief Calculate Hydrodynamic Fluxes using the Riemann solver

void Hydro::CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                            AthenaArray<Real> &bcc, const int order) {
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;

  // b,bcc are passed as fn parameters becausse clients may want to pass different bcc1,
  // b1, b2, etc., but the remaining members of the Field class are accessed directly via
  // pointers because they are unique. NOTE: b, bcc are nullptrs if no MHD.
#if MAGNETIC_FIELDS_ENABLED
  // used only to pass to (up-to) 2x RiemannSolver() calls per dimension:
  // x1:
  AthenaArray<Real> &b1 = b.x1f, &w_x1f = pmb->pfield->wght.x1f,
                  &e3x1 = pmb->pfield->e3_x1f, &e2x1 = pmb->pfield->e2_x1f;
  // x2:
  AthenaArray<Real> &b2 = b.x2f, &w_x2f = pmb->pfield->wght.x2f,
                  &e1x2 = pmb->pfield->e1_x2f, &e3x2 = pmb->pfield->e3_x2f;
  // x3:
  AthenaArray<Real> &b3 = b.x3f, &w_x3f = pmb->pfield->wght.x3f,
                  &e1x3 = pmb->pfield->e1_x3f, &e2x3 = pmb->pfield->e2_x3f;
#endif
  AthenaArray<Real> &flux_fc = scr1_nkji_;
  AthenaArray<Real> &laplacian_all_fc = scr2_nkji_;

  //--------------------------------------------------------------------------------------
  // i-direction

  AthenaArray<Real> &x1flux = flux[X1DIR];
  // set the loop limits
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED || order == 4) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl = js-1, ju = je+1, kl = ks, ku = ke;
      else // 3D
        jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if (order == 1) {
        pmb->precon->DonorCellX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
        else
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, w, bcc, wl_, wr_);
      }

#if EOS_SCALAR_INPUT_ENABLED
      if (NSCALARS > 0) {
        AthenaArray<Real> &r = pmb->pscalars->r;
        
        if (order == 1) {
          pmb->precon->DonorCellX1(
              k, j, is-1, ie+1, r, rl_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX1(
              k, j, is-1, ie+1, r, rl_, rr_);
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in Hydro::CalculateFluxes" << std::endl
              << "EOS scalar input currently supports reconstruction order <= 2."
              << std::endl;
          ATHENA_ERROR(msg);
        }
      }
#endif

      pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
#if EOS_SCALAR_INPUT_ENABLED
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_, &rl_, &rr_);
#else
      RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, x1flux, dxw_);
#endif
#else  // MHD:
      // x1flux(IBY) = (v1*b2 - v2*b1) = -EMFZ
      // x1flux(IBZ) = (v1*b3 - v3*b1) =  EMFY
      RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, x1flux, e3x1, e2x1, w_x1f, dxw_);
#endif

      if (order == 4) {
        for (int n=0; n<NWAVE; n++) {
          for (int i=is; i<=ie+1; i++) {
            wl3d_(n,k,j,i) = wl_(n,i);
            wr3d_(n,k,j,i) = wr_(n,i);
          }
        }
      }
    }
  }

  if (order == 4) {
    // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so this should factor out
    // TODO(felker): also, this may need to be dx1v, since Laplacian is cell-centered
    Real h = pmb->pcoord->dx1f(is);  // pco->dx1f(i); inside loop
    Real C = (h*h)/24.0;

    // construct Laplacian from x1flux
    pmb->pcoord->LaplacianX1All(x1flux, laplacian_all_fc, 0, NHYDRO-1,
                                kl, ku, jl, ju, is, ie+1);

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // Compute Laplacian of primitive Riemann states on x1 faces
        for (int n=0; n<NWAVE; ++n) {
          pmb->pcoord->LaplacianX1(wl3d_, laplacian_l_fc_, n, k, j, is, ie+1);
          pmb->pcoord->LaplacianX1(wr3d_, laplacian_r_fc_, n, k, j, is, ie+1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
            wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
          pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
        }

        // Compute x1 interface fluxes from face-centered primitive variables
        // TODO(felker): check that e3x1,e2x1 arguments added in late 2017 work here
        pmb->pcoord->CenterWidth1(k, j, is, ie+1, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
        RiemannSolver(k, j, is, ie+1, IVX, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
        RiemannSolver(k, j, is, ie+1, IVX, b1, wl_, wr_, flux_fc, e3x1, e2x1,
                      w_x1f, dxw_);
#endif
        // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie+1; i++) {
            x1flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
            // TODO(felker): replace this loop-based deep copy with memcpy, or alternative
            if (n == IDN && NSCALARS > 0) {
              pmb->pscalars->mass_flux_fc[X1DIR](k,j,i) = flux_fc(n,k,j,i);
            }
          }
        }
      }
    }
  } // end if (order == 4)
  //------------------------------------------------------------------------------
  // end x1 fourth-order hydro

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmb->pmy_mesh->f2) {
    AthenaArray<Real> &x2flux = flux[X2DIR];
    // set the loop limits
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      if (pmb->block_size.nx3 == 1) // 2D
        kl = ks, ku = ke;
      else // 3D
        kl = ks-1, ku = ke+1;
    }

    for (int k=kl; k<=ku; ++k) {
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX2(k, js-1, il, iu, w, bcc, wl_, wr_);
        else
          pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, w, bcc, wl_, wr_);
      }
#if EOS_SCALAR_INPUT_ENABLED
      AthenaArray<Real> &r = pmb->pscalars->r;
      
      if (order == 1) {
        pmb->precon->DonorCellX2(k, js-1, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, r, rl_, rr_);
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in Hydro::CalculateFluxes" << std::endl
            << "EOS scalar input currently supports reconstruction order <= 2."
            << std::endl;
        ATHENA_ERROR(msg);
      }
#endif

      for (int j=js; j<=je+1; ++j) {
        // reconstruct L/R states at j
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
            pmb->precon->PiecewiseParabolicFastX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
            pmb->precon->WENOZX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
            pmb->precon->WENOMZX2(k, j, il, iu, w, bcc, wlb_, wr_);
          else
            pmb->precon->PiecewiseParabolicX2(k, j, il, iu, w, bcc, wlb_, wr_);
        }

        // Capture the actual PLM output and limiter scratch arrays before any
        // scalar reconstruction reuses Reconstruction's scratch storage.  This
        // distinguishes a bad hydro slope from later buffer corruption.
        if (pmb->gid == 0 && k == ks && order == 2 && j >= js && j <= js+1) {
          const int side = j-js;
          const int vars[2] = {IDN, IPR};
          for (int nv=0; nv<2; ++nv) {
            const int n = vars[nv];
            polar_debug_.face_plm_after_hydro[side][nv][0] = wlb_(n,is);
            polar_debug_.face_plm_after_hydro[side][nv][1] = wr_(n,is);
            polar_debug_.face_plm_scratch[side][nv][0] = pmb->precon->scr1_ni_(n,is);
            polar_debug_.face_plm_scratch[side][nv][1] = pmb->precon->scr2_ni_(n,is);
            polar_debug_.face_plm_scratch[side][nv][2] = pmb->precon->scr3_ni_(n,is);
            polar_debug_.face_plm_scratch[side][nv][3] = pmb->precon->scr4_ni_(n,is);
          }
        }
#if EOS_SCALAR_INPUT_ENABLED
        // scalar: current row
        if (order == 1) {
          pmb->precon->DonorCellX2(k, j, il, iu, r, rlb_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, j, il, iu, r, rlb_, rr_);
        }
#endif

        if (pmb->gid == 0 && k == ks && order == 2 && j >= js && j <= js+1) {
          const int side = j-js;
          const int vars[2] = {IDN, IPR};
          for (int nv=0; nv<2; ++nv) {
            const int n = vars[nv];
            polar_debug_.face_plm_after_scalars[side][nv][0] = wlb_(n,is);
            polar_debug_.face_plm_after_scalars[side][nv][1] = wr_(n,is);
          }
        }

        pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);

        // Preserve the exact untransformed PLM inputs and reconstruct the scalar
        // limiter arithmetic for the two cells adjacent to the known failing face.
        // Nothing is printed here; AddFluxDivergence reports this buffer only if the
        // corresponding active-cell density update is non-positive.
        if (pmb->gid == 0 && k == ks && j == js+1 && order == 2) {
          Coordinates *pco = pmb->pcoord;
          polar_debug_.face_order = order;
          polar_debug_.face_characteristic_projection =
              pmb->precon->characteristic_projection_;
          polar_debug_.face_minmod = pmb->precon->minmod_;
          for (int n=0; n<NWAVE; ++n) {
            polar_debug_.face_wl_pre[n] = wl_(n,is);
            polar_debug_.face_wr_pre[n] = wr_(n,is);
          }

          const int vars[2] = {IDN, IPR};
          const int cells[2] = {js, js+1};
          for (int side=0; side<2; ++side) {
            const int jc = cells[side];
            const Real cf = pco->dx2v(jc)/(pco->x2f(jc+1)-pco->x2v(jc));
            const Real cb = pco->dx2v(jc-1)/(pco->x2v(jc)-pco->x2f(jc));
            const Real dxF = pco->dx2f(jc)/pco->dx2v(jc);
            const Real dxB = pco->dx2f(jc)/pco->dx2v(jc-1);
            const Real dxp = (pco->x2f(jc+1)-pco->x2v(jc))/pco->dx2f(jc);
            const Real dxm = (pco->x2v(jc)-pco->x2f(jc))/pco->dx2f(jc);
            polar_debug_.face_stencil_geom[side][0] = cf;
            polar_debug_.face_stencil_geom[side][1] = cb;
            polar_debug_.face_stencil_geom[side][2] = dxF;
            polar_debug_.face_stencil_geom[side][3] = dxB;
            polar_debug_.face_stencil_geom[side][4] = dxp;
            polar_debug_.face_stencil_geom[side][5] = dxm;

            for (int nv=0; nv<2; ++nv) {
              const int n = vars[nv];
              const Real qm = w(n,k,jc-1,is);
              const Real qc = w(n,k,jc,is);
              const Real qp = w(n,k,jc+1,is);
              const Real dwl = qc-qm;
              const Real dwr = qp-qc;
              const Real dqB = dwl*dxB;
              const Real dqF = dwr*dxF;
              Real dwm = 0.0;
              if (!pmb->precon->characteristic_projection_) {
                if (pmb->precon->minmod_) {
                  if (dwl*dwr > 0.0) {
                    const Real dwlw = dwl*dxB;
                    const Real dwrw = dwr*dxF;
                    dwm = (dwlw >= 0.0) ? std::min(dwlw,dwrw)
                                         : std::max(dwlw,dwrw);
                  }
                } else {
                  const Real dq2 = dqF*dqB;
                  dwm = dq2*(cf*dqB+cb*dqF)
                        /(SQR(dqB)+SQR(dqF)+dq2*(cf+cb-2.0));
                  if (dq2 <= 0.0) dwm = 0.0;
                }
              }
              polar_debug_.face_stencil_q[side][nv][0] = qm;
              polar_debug_.face_stencil_q[side][nv][1] = qc;
              polar_debug_.face_stencil_q[side][nv][2] = qp;
              polar_debug_.face_stencil_calc[side][nv][0] = dwl;
              polar_debug_.face_stencil_calc[side][nv][1] = dwr;
              polar_debug_.face_stencil_calc[side][nv][2] = dqB;
              polar_debug_.face_stencil_calc[side][nv][3] = dqF;
              polar_debug_.face_stencil_calc[side][nv][4] = dwm;
              polar_debug_.face_stencil_calc[side][nv][5] = qc-dxm*dwm;
              polar_debug_.face_stencil_calc[side][nv][6] = qc+dxp*dwm;
            }
          }
        }
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
#if EOS_SCALAR_INPUT_ENABLED
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_, &rl_, &rr_);
#else
        RiemannSolver(k, j, il, iu, IVY, wl_, wr_, x2flux, dxw_);
#endif
#else  // MHD:
        // flx(IBY) = (v2*b3 - v3*b2) = -EMFX
        // flx(IBZ) = (v2*b1 - v1*b2) =  EMFZ
        RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, x2flux, e1x2, e3x2, w_x2f, dxw_);
#endif

        // Save the exact, locally transformed HLLC inputs at the upper x2 face of the
        // first active polar cell.  RiemannSolver() transforms wl_/wr_ in place in GR.
        if (pmb->gid == 0 && k == ks && j == js+1) {
          polar_debug_.face_valid = true;
          polar_debug_.face_cycle = pmb->pmy_mesh->ncycle;
          polar_debug_.face_time = pmb->pmy_mesh->time;
          for (int n=0; n<NWAVE; ++n) {
            polar_debug_.face_wl[n] = wl_(n,is);
            polar_debug_.face_wr[n] = wr_(n,is);
          }

          const Real gamma_prime = pmb->peos->GetGamma()
                                   /(pmb->peos->GetGamma() - 1.0);
          const Real ul0 = std::sqrt(1.0 + SQR(wl_(IVX,is)) + SQR(wl_(IVY,is))
                                     + SQR(wl_(IVZ,is)));
          const Real ur0 = std::sqrt(1.0 + SQR(wr_(IVX,is)) + SQR(wr_(IVY,is))
                                     + SQR(wr_(IVZ,is)));
          const Real wgas_l = wl_(IDN,is) + gamma_prime*wl_(IPR,is);
          const Real wgas_r = wr_(IDN,is) + gamma_prime*wr_(IPR,is);
          polar_debug_.hllc_wgas_l = wgas_l;
          polar_debug_.hllc_wgas_r = wgas_r;
          polar_debug_.hllc_cs2_l = pmb->peos->GetGamma()*wl_(IPR,is)/wgas_l;
          polar_debug_.hllc_cs2_r = pmb->peos->GetGamma()*wr_(IPR,is)/wgas_r;
          polar_debug_.hllc_cons_d_l = wl_(IDN,is)*ul0;
          polar_debug_.hllc_cons_d_r = wr_(IDN,is)*ur0;
          pmb->peos->SoundSpeedsSR(wgas_l, wl_(IPR,is), wl_(IVY,is)/ul0, SQR(ul0),
                                   &polar_debug_.lambda_p_l,
                                   &polar_debug_.lambda_m_l);
          pmb->peos->SoundSpeedsSR(wgas_r, wr_(IPR,is), wr_(IVY,is)/ur0, SQR(ur0),
                                   &polar_debug_.lambda_p_r,
                                   &polar_debug_.lambda_m_r);

          const Real lambda_l = std::min(polar_debug_.lambda_m_l,
                                         polar_debug_.lambda_m_r);
          const Real lambda_r = std::max(polar_debug_.lambda_p_l,
                                         polar_debug_.lambda_p_r);
          polar_debug_.hllc_lambda_l = lambda_l;
          polar_debug_.hllc_lambda_r = lambda_r;
          Real cons_l[NWAVE], cons_r[NWAVE], flux_l[NWAVE], flux_r[NWAVE];
          const int ivx = IVY;
          const int ivy = IVX + ((ivx-IVX)+1)%3;
          const int ivz = IVX + ((ivx-IVX)+2)%3;
          for (int n=0; n<NWAVE; ++n) {
            cons_l[n] = cons_r[n] = flux_l[n] = flux_r[n] = 0.0;
          }
          cons_l[IDN] = wl_(IDN,is)*ul0;
          cons_l[IEN] = wgas_l*SQR(ul0)-wl_(IPR,is);
          cons_l[ivx] = wgas_l*wl_(ivx,is)*ul0;
          cons_l[ivy] = wgas_l*wl_(ivy,is)*ul0;
          cons_l[ivz] = wgas_l*wl_(ivz,is)*ul0;
          flux_l[IDN] = wl_(IDN,is)*wl_(ivx,is);
          flux_l[IEN] = wgas_l*ul0*wl_(ivx,is);
          flux_l[ivx] = wgas_l*SQR(wl_(ivx,is))+wl_(IPR,is);
          flux_l[ivy] = wgas_l*wl_(ivy,is)*wl_(ivx,is);
          flux_l[ivz] = wgas_l*wl_(ivz,is)*wl_(ivx,is);
          cons_r[IDN] = wr_(IDN,is)*ur0;
          cons_r[IEN] = wgas_r*SQR(ur0)-wr_(IPR,is);
          cons_r[ivx] = wgas_r*wr_(ivx,is)*ur0;
          cons_r[ivy] = wgas_r*wr_(ivy,is)*ur0;
          cons_r[ivz] = wgas_r*wr_(ivz,is)*ur0;
          flux_r[IDN] = wr_(IDN,is)*wr_(ivx,is);
          flux_r[IEN] = wgas_r*ur0*wr_(ivx,is);
          flux_r[ivx] = wgas_r*SQR(wr_(ivx,is))+wr_(IPR,is);
          flux_r[ivy] = wgas_r*wr_(ivy,is)*wr_(ivx,is);
          flux_r[ivz] = wgas_r*wr_(ivz,is)*wr_(ivx,is);
          const Real inv_dlambda = 1.0/(lambda_r-lambda_l);
          Real cons_hll[NWAVE], flux_hll[NWAVE];
          for (int n=0; n<NWAVE; ++n) {
            cons_hll[n] = (lambda_r*cons_r[n]-lambda_l*cons_l[n]
                           +flux_l[n]-flux_r[n])*inv_dlambda;
            flux_hll[n] = (lambda_r*flux_l[n]-lambda_l*flux_r[n]
                           +lambda_l*lambda_r*(cons_r[n]-cons_l[n]))*inv_dlambda;
          }
          const Real contact_b = -(cons_hll[IEN]+flux_hll[ivx]);
          const Real contact_disc = SQR(contact_b)
                                    -4.0*flux_hll[IEN]*cons_hll[ivx];
          polar_debug_.hllc_contact_discriminant = contact_disc;
          Real lambda_star;
          if (std::abs(flux_hll[IEN]) > TINY_NUMBER) {
            lambda_star = -2.0*cons_hll[ivx]
                          /(contact_b-std::sqrt(contact_disc));
          } else {
            lambda_star = -cons_hll[ivx]/contact_b;
          }
          polar_debug_.hllc_lambda_star = lambda_star;
          polar_debug_.hllc_pgas_star =
              -flux_hll[IEN]*lambda_star+flux_hll[ivx];
          polar_debug_.hllc_cons_d_lstar =
              cons_l[IDN]*(lambda_l-wl_(ivx,is)/ul0)/(lambda_l-lambda_star);
          polar_debug_.hllc_cons_d_rstar =
              cons_r[IDN]*(lambda_r-wr_(ivx,is)/ur0)/(lambda_r-lambda_star);
        }

        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
#if EOS_SCALAR_INPUT_ENABLED
        rl_.SwapAthenaArray(rlb_);
#endif
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx2v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx2f(js);  // pco->dx2f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x2flux
      pmb->pcoord->LaplacianX2All(x2flux, laplacian_all_fc, 0, NHYDRO-1,
                                  kl, ku, js, je+1, il, iu);

      // Approximate x2 face-centered primitive Riemann states
      for (int k=kl; k<=ku; ++k) {
        for (int j=js; j<=je+1; ++j) {
          // Compute Laplacian of primitive Riemann states on x2 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX2(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX2(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
          }

          // Compute x2 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e1x2,e3x2 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth2(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVY, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVY, b2, wl_, wr_, flux_fc, e1x2, e3x2,
                        w_x2f, dxw_);
#endif

          // Apply Laplacian of second-order accurate face-averaged flux on x1 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x2flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
              if (n == IDN && NSCALARS > 0) {
                pmb->pscalars->mass_flux_fc[X2DIR](k,j,i) = flux_fc(n,k,j,i);
              }
            }
          }
        }
      }
    } // end if (order == 4)
  }

  //--------------------------------------------------------------------------------------
  // k-direction

  if (pmb->pmy_mesh->f3) {
    AthenaArray<Real> &x3flux = flux[X3DIR];
    // set the loop limits
    il = is, iu = ie, jl = js, ju = je;
    if (MAGNETIC_FIELDS_ENABLED || order == 4) {
      il = is-1, iu = ie+1, jl = js-1, ju = je+1;
    }

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if (order == 1) {
        pmb->precon->DonorCellX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      } else {
        if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
          pmb->precon->PiecewiseParabolicFastX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
          pmb->precon->WENOZX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
          pmb->precon->WENOMZX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
        else
          pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, w, bcc, wl_, wr_);
      }
#if EOS_SCALAR_INPUT_ENABLED
      AthenaArray<Real> &r = pmb->pscalars->r;

      if (order == 1) {
        pmb->precon->DonorCellX3(
            ks-1, j, il, iu, r, rl_, rr_);
      } else if (order == 2) {
        pmb->precon->PiecewiseLinearX3(
            ks-1, j, il, iu, r, rl_, rr_);
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in Hydro::CalculateFluxes" << std::endl
            << "EOS scalar input currently supports reconstruction order <= 2."
            << std::endl;
        ATHENA_ERROR(msg);
      }
#endif
      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if (order == 1) {
          pmb->precon->DonorCellX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(k, j, il, iu, w, bcc, wlb_, wr_);
        } else {
          if (pmb->precon->rec3m_ == REC3METHOD::PPMF)
            pmb->precon->PiecewiseParabolicFastX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOZ)
            pmb->precon->WENOZX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else if (pmb->precon->rec3m_ == REC3METHOD::WENOMZ)
            pmb->precon->WENOMZX3(k, j, il, iu, w, bcc, wlb_, wr_);
          else
            pmb->precon->PiecewiseParabolicX3(k, j, il, iu, w, bcc, wlb_, wr_);
        }
#if EOS_SCALAR_INPUT_ENABLED
        if (order == 1) {
          pmb->precon->DonorCellX3(
              k, j, il, iu, r, rlb_, rr_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(
              k, j, il, iu, r, rlb_, rr_);
        }
#endif
        
        pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
#if EOS_SCALAR_INPUT_ENABLED
        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_, &rl_, &rr_);
#else
        RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, x3flux, dxw_);
#endif
#else  // MHD:
        // flx(IBY) = (v3*b1 - v1*b3) = -EMFY
        // flx(IBZ) = (v3*b2 - v2*b3) =  EMFX
        RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, x3flux, e2x3, e1x3, w_x3f, dxw_);
#endif
        if (order == 4) {
          for (int n=0; n<NWAVE; n++) {
            for (int i=il; i<=iu; i++) {
              wl3d_(n,k,j,i) = wl_(n,i);
              wr3d_(n,k,j,i) = wr_(n,i);
            }
          }
        }

        // swap the arrays for the next step
        wl_.SwapAthenaArray(wlb_);
#if EOS_SCALAR_INPUT_ENABLED
        rl_.SwapAthenaArray(rlb_);
#endif
      }
    }
    if (order == 4) {
      // TODO(felker): assuming uniform mesh with dx1f=dx2f=dx3f, so factor this out
      // TODO(felker): also, this may need to be dx3v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx3f(ks);  // pco->dx3f(j); inside loop
      Real C = (h*h)/24.0;

      // construct Laplacian from x3flux
      pmb->pcoord->LaplacianX3All(x3flux, laplacian_all_fc, 0, NHYDRO-1,
                                  ks, ke+1, jl, ju, il, iu);

      // Approximate x3 face-centered primitive Riemann states
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=jl; j<=ju; ++j) {
          // Compute Laplacian of primitive Riemann states on x3 faces
          for (int n=0; n<NWAVE; ++n) {
            pmb->pcoord->LaplacianX3(wl3d_, laplacian_l_fc_, n, k, j, il, iu);
            pmb->pcoord->LaplacianX3(wr3d_, laplacian_r_fc_, n, k, j, il, iu);
#pragma omp simd
            for (int i=il; i<=iu; ++i) {
              wl_(n,i) = wl3d_(n,k,j,i) - C*laplacian_l_fc_(i);
              wr_(n,i) = wr3d_(n,k,j,i) - C*laplacian_r_fc_(i);
            }
          }
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            pmb->peos->ApplyPrimitiveFloors(wl_, k, j, i);
            pmb->peos->ApplyPrimitiveFloors(wr_, k, j, i);
          }

          // Compute x3 interface fluxes from face-centered primitive variables
          // TODO(felker): check that e2x3,e1x3 arguments added in late 2017 work here
          pmb->pcoord->CenterWidth3(k, j, il, iu, dxw_);
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
          RiemannSolver(k, j, il, iu, IVZ, wl_, wr_, flux_fc, dxw_);
#else  // MHD:
          RiemannSolver(k, j, il, iu, IVZ, b3, wl_, wr_, flux_fc, e2x3, e1x3,
                        w_x3f, dxw_);
#endif
          // Apply Laplacian of second-order accurate face-averaged flux on x3 faces
          for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
            for (int i=il; i<=iu; i++) {
              x3flux(n,k,j,i) = flux_fc(n,k,j,i) + C*laplacian_all_fc(n,k,j,i);
              if (n == IDN && NSCALARS > 0) {
                pmb->pscalars->mass_flux_fc[X3DIR](k,j,i) = flux_fc(n,k,j,i);
              }
            }
          }
        }
      }
    } // end if (order == 4)
  }

  if (!STS_ENABLED)
    AddDiffusionFluxes();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalculateFluxes_STS
//! \brief Calculate Hydrodynamic Diffusion Fluxes for STS

void Hydro::CalculateFluxes_STS() {
  AddDiffusionFluxes();
}

void Hydro::AddDiffusionFluxes() {
  Field *pf = pmy_block->pfield;
  // add diffusion fluxes
  if (hdif.hydro_diffusion_defined) {
    if (hdif.nu_iso > 0.0 || hdif.nu_aniso > 0.0)
      hdif.AddDiffusionFlux(hdif.visflx,flux);
    if (NON_BAROTROPIC_EOS) {
      if (hdif.kappa_iso > 0.0 || hdif.kappa_aniso > 0.0)
        hdif.AddDiffusionEnergyFlux(hdif.cndflx,flux);
    }
  }
  if (MAGNETIC_FIELDS_ENABLED && NON_BAROTROPIC_EOS) {
    if (pf->fdif.field_diffusion_defined)
      pf->fdif.AddPoyntingFlux(pf->fdif.pflux);
  }
  return;
}

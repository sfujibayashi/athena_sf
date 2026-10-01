//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file add_flux_divergence.cpp
//! \brief Computes divergence of the Hydro fluxes and adds that to a conserved variable
//! register

// C headers
#include <cstdlib>    // std::abort

// C++ headers
#include <algorithm>  // std::binary_search
#include <vector>     // std::vector
#include <cmath>      // std::abs, std::sqrt
#include <iomanip>    // std::setprecision
#include <limits>     // std::numeric_limits

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "hydro.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AddFluxDivergence
//! \brief Adds flux divergence to weighted average of conservative variables from
//! previous step(s) of time integrator algorithm

// TODO(felker): consider combining with PassiveScalars implementation + (see 57cfe28b)
// (may rename to AddPhysicalFluxDivergence or AddQuantityFluxDivergence to explicitly
// distinguish from CoordTerms)
// (may rename to AddHydroFluxDivergence and AddScalarsFluxDivergence, if
// the implementations remain completely independent / no inheritance is
// used)
void Hydro::AddFluxDivergence(const Real wght, AthenaArray<Real> &u_out) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = flux[X1DIR];
  AthenaArray<Real> &x2flux = flux[X2DIR];
  AthenaArray<Real> &x3flux = flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, is, ie+1, x1area);
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = (x1area(i+1)*x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j  , is, ie, x2area   );
        pmb->pcoord->Face2Area(k, j+1, is, ie, x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  , j, is, ie, x3area   );
        pmb->pcoord->Face3Area(k+1, j, is, ie, x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k, j, is, ie, vol);

      // Target only the known first active polar cell.  Keeping this outside the
      // conserved-variable update loop avoids perturbing the arithmetic of the update.
      if (pmb->gid == 0 && k == ks && j == js) {
        const int i = is;
        const Real d_old = u_out(IDN,k,j,i);
        const Real df1 = x1area(i+1)*x1flux(IDN,k,j,i+1)
                         - x1area(i)*x1flux(IDN,k,j,i);
        const Real df2 = x2area_p1(i)*x2flux(IDN,k,j+1,i)
                         - x2area(i)*x2flux(IDN,k,j,i);
        const Real dD1 = -wght*df1/vol(i);
        const Real dD2 = -wght*df2/vol(i);
        const Real d_new = d_old + dD1 + dD2;

        if (d_new <= 0.0) {
          Coordinates *pco = pmb->pcoord;
          Mesh *pm = pmb->pmy_mesh;
          const Real inf = std::numeric_limits<Real>::infinity();
          const Real dt_pos_x1 = (df1 > 0.0) ? d_old*vol(i)/df1 : inf;
          const Real dt_pos_x2 = (df2 > 0.0) ? d_old*vol(i)/df2 : inf;
          const Real df_total = df1 + df2;
          const Real dt_pos_total = (df_total > 0.0) ? d_old*vol(i)/df_total : inf;

          // Current cell metric and coordinate light speeds.
          pco->CellMetric(k, j, i, i, g_, gi_);
          const Real g00 = g_(I00,i);
          const Real g02 = g_(I02,i);
          const Real g22 = g_(I22,i);
          const Real gi00 = gi_(I00,i);
          const Real gi02 = gi_(I02,i);
          const Real gi22 = gi_(I22,i);
          const Real alpha = std::sqrt(-1.0/gi00);
          const Real beta2 = -gi02/gi00;
          const Real speed1_cell =
              -(std::sqrt(SQR(gi_(I01,i)) - gi00*gi_(I11,i))
                + std::abs(gi_(I01,i)))/gi00;
          const Real speed2_cell =
              -(std::sqrt(SQR(gi02) - gi00*gi22) + std::abs(gi02))/gi00;

          // Face coordinate light speeds, paired with the exact A/V factors used above.
          pco->Face1Metric(k, j, i, i, g_, gi_);
          const Real speed1_l =
              -(std::sqrt(SQR(gi_(I01,i)) - gi_(I00,i)*gi_(I11,i))
                + std::abs(gi_(I01,i)))/gi_(I00,i);
          pco->Face1Metric(k, j, i+1, i+1, g_, gi_);
          const Real speed1_r =
              -(std::sqrt(SQR(gi_(I01,i+1)) - gi_(I00,i+1)*gi_(I11,i+1))
                + std::abs(gi_(I01,i+1)))/gi_(I00,i+1);
          Real speed2_l = 0.0;
          if (x2area(i) != 0.0) {
            pco->Face2Metric(k, j, i, i, g_, gi_);
            speed2_l = -(std::sqrt(SQR(gi_(I02,i)) - gi_(I00,i)*gi_(I22,i))
                         + std::abs(gi_(I02,i)))/gi_(I00,i);
          }
          pco->Face2Metric(k, j+1, i, i, g_, gi_);
          const Real speed2_r =
              -(std::sqrt(SQR(gi_(I02,i)) - gi_(I00,i)*gi_(I22,i))
                + std::abs(gi_(I02,i)))/gi_(I00,i);

          const Real rate1 = std::max(x1area(i)*speed1_l,
                                      x1area(i+1)*speed1_r)/vol(i);
          const Real rate2 = std::max(x2area(i)*speed2_l,
                                      x2area_p1(i)*speed2_r)/vol(i);
          const Real fv_dt1_raw = 1.0/rate1;
          const Real fv_dt2_raw = 1.0/rate2;
          const Real fv_dt_multid_raw = 1.0/(rate1 + rate2);
          const Real cfl = pm->cfl_number;

          std::cout << std::setprecision(16)
                    << "\nPOLAR_D_NEGATIVE_BEGIN\n"
                    << "cycle=" << pm->ncycle
                    << " time=" << pm->time
                    << " dt=" << pm->dt
                    << " gid=" << pmb->gid
                    << " i=" << i << " j=" << j << " k=" << k
                    << " wght=" << wght
                    << " wght_over_dt=" << wght/pm->dt << "\n"
                    << "D_old=" << d_old
                    << " D_new=" << d_new
                    << " dD_x1=" << dD1
                    << " dD_x2=" << dD2 << "\n"
                    << "F1L_D=" << x1flux(IDN,k,j,i)
                    << " F1R_D=" << x1flux(IDN,k,j,i+1)
                    << " F2L_D=" << x2flux(IDN,k,j,i)
                    << " F2R_D=" << x2flux(IDN,k,j+1,i) << "\n"
                    << "df1=" << df1 << " df2=" << df2
                    << " df_total=" << df_total << "\n"
                    << "A1L=" << x1area(i) << " A1R=" << x1area(i+1)
                    << " A2L=" << x2area(i) << " A2R=" << x2area_p1(i)
                    << " V=" << vol(i)
                    << " V_over_A2R=" << vol(i)/x2area_p1(i) << "\n"
                    << "rho=" << w(IDN,k,j,i)
                    << " p=" << w(IPR,k,j,i)
                    << " u1=" << w(IVX,k,j,i)
                    << " u2=" << w(IVY,k,j,i)
                    << " u3=" << w(IVZ,k,j,i) << "\n"
                    << "metric_current: g00=" << g00
                    << " g02=" << g02 << " g22=" << g22
                    << " gi00=" << gi00 << " gi02=" << gi02
                    << " gi22=" << gi22
                    << " alpha=" << alpha << " beta2=" << beta2 << "\n"
                    << "light_speeds: cell_x1=" << speed1_cell
                    << " cell_x2=" << speed2_cell
                    << " face_x1L=" << speed1_l
                    << " face_x1R=" << speed1_r
                    << " face_x2L=" << speed2_l
                    << " face_x2R=" << speed2_r << "\n"
                    << "fv_rates: rate1=" << rate1 << " rate2=" << rate2
                    << " sum=" << rate1+rate2 << "\n"
                    << "fv_dt_raw: x1=" << fv_dt1_raw
                    << " x2=" << fv_dt2_raw
                    << " multid=" << fv_dt_multid_raw << "\n"
                    << "fv_dt_cfl: x1=" << cfl*fv_dt1_raw
                    << " x2=" << cfl*fv_dt2_raw
                    << " multid=" << cfl*fv_dt_multid_raw << "\n"
                    << "dt_pos: x1=" << dt_pos_x1
                    << " x2=" << dt_pos_x2
                    << " total=" << dt_pos_total << "\n";

          if (polar_debug_.cfl_valid) {
            std::cout << "saved_NewBlockTimeStep: computed_cycle="
                      << polar_debug_.cfl_cycle
                      << " computed_time=" << polar_debug_.cfl_time
                      << " r=" << polar_debug_.r
                      << " theta=" << polar_debug_.theta
                      << " dtheta=" << polar_debug_.dtheta << "\n"
                      << "saved_widths: CenterWidth1=" << polar_debug_.center_width1
                      << " CenterWidth2=" << polar_debug_.center_width2 << "\n"
                      << "saved_metric: g00=" << polar_debug_.g00
                      << " g02=" << polar_debug_.g02
                      << " g22=" << polar_debug_.g22
                      << " gi00=" << polar_debug_.gi00
                      << " gi02=" << polar_debug_.gi02
                      << " gi22=" << polar_debug_.gi22 << "\n"
                      << "saved_speeds: speed1=" << polar_debug_.speed1
                      << " speed2=" << polar_debug_.speed2 << "\n"
                      << "saved_dt_raw: target_x1=" << polar_debug_.dt1_raw
                      << " target_x2=" << polar_debug_.dt2_raw
                      << " block_min=" << polar_debug_.block_dt_raw << "\n"
                      << "saved_dt_cfl: target_x1=" << cfl*polar_debug_.dt1_raw
                      << " target_x2=" << cfl*polar_debug_.dt2_raw
                      << " block_min=" << polar_debug_.block_dt_scaled << "\n"
                      << "saved_fv_dt_raw: x1=" << polar_debug_.fv_dt1_raw
                      << " x2=" << polar_debug_.fv_dt2_raw
                      << " saved_fv_multid_cfl="
                      << polar_debug_.fv_dt_multid_scaled << "\n"
                      << "saved_primitive: rho=" << polar_debug_.rho
                      << " p=" << polar_debug_.press
                      << " u1=" << polar_debug_.u1
                      << " u2=" << polar_debug_.u2
                      << " u3=" << polar_debug_.u3 << "\n"
                      << "saved_geometry: A1L=" << polar_debug_.area1_l
                      << " A1R=" << polar_debug_.area1_r
                      << " A2L=" << polar_debug_.area2_l
                      << " A2R=" << polar_debug_.area2_r
                      << " V=" << polar_debug_.volume << "\n";
          }

          if (polar_debug_.face_valid) {
            std::cout << "x2_face_state: cycle=" << polar_debug_.face_cycle
                      << " time=" << polar_debug_.face_time
                      << " order=" << polar_debug_.face_order
                      << " characteristic_projection="
                      << polar_debug_.face_characteristic_projection
                      << " minmod=" << polar_debug_.face_minmod << "\n"
                      << "WL_global_pre: rho=" << polar_debug_.face_wl_pre[IDN]
                      << " p=" << polar_debug_.face_wl_pre[IPR]
                      << " u1=" << polar_debug_.face_wl_pre[IVX]
                      << " u2=" << polar_debug_.face_wl_pre[IVY]
                      << " u3=" << polar_debug_.face_wl_pre[IVZ] << "\n"
                      << "WR_global_pre: rho=" << polar_debug_.face_wr_pre[IDN]
                      << " p=" << polar_debug_.face_wr_pre[IPR]
                      << " u1=" << polar_debug_.face_wr_pre[IVX]
                      << " u2=" << polar_debug_.face_wr_pre[IVY]
                      << " u3=" << polar_debug_.face_wr_pre[IVZ] << "\n"
                      << "WL_local: rho=" << polar_debug_.face_wl[IDN]
                      << " p=" << polar_debug_.face_wl[IPR]
                      << " u1=" << polar_debug_.face_wl[IVX]
                      << " u2=" << polar_debug_.face_wl[IVY]
                      << " u3=" << polar_debug_.face_wl[IVZ] << "\n"
                      << "WR_local: rho=" << polar_debug_.face_wr[IDN]
                      << " p=" << polar_debug_.face_wr[IPR]
                      << " u1=" << polar_debug_.face_wr[IVX]
                      << " u2=" << polar_debug_.face_wr[IVY]
                      << " u3=" << polar_debug_.face_wr[IVZ] << "\n"
                      << "x2_lambdas_local: lp_l=" << polar_debug_.lambda_p_l
                      << " lm_l=" << polar_debug_.lambda_m_l
                      << " lp_r=" << polar_debug_.lambda_p_r
                      << " lm_r=" << polar_debug_.lambda_m_r << "\n"
                      << "hllc_thermo: wgas_l=" << polar_debug_.hllc_wgas_l
                      << " wgas_r=" << polar_debug_.hllc_wgas_r
                      << " cs2_l=" << polar_debug_.hllc_cs2_l
                      << " cs2_r=" << polar_debug_.hllc_cs2_r
                      << " D_l=" << polar_debug_.hllc_cons_d_l
                      << " D_r=" << polar_debug_.hllc_cons_d_r << "\n"
                      << "hllc_star: lambda_l=" << polar_debug_.hllc_lambda_l
                      << " lambda_r=" << polar_debug_.hllc_lambda_r
                      << " lambda_star=" << polar_debug_.hllc_lambda_star
                      << " pgas_star=" << polar_debug_.hllc_pgas_star
                      << " discriminant=" << polar_debug_.hllc_contact_discriminant
                      << " D_lstar=" << polar_debug_.hllc_cons_d_lstar
                      << " D_rstar=" << polar_debug_.hllc_cons_d_rstar << "\n";
            const char *side_name[2] = {"left_cell_j2", "right_cell_j3"};
            const char *var_name[2] = {"rho", "p"};
            for (int side=0; side<2; ++side) {
              std::cout << "plm_geom " << side_name[side]
                        << ": cf=" << polar_debug_.face_stencil_geom[side][0]
                        << " cb=" << polar_debug_.face_stencil_geom[side][1]
                        << " dxF=" << polar_debug_.face_stencil_geom[side][2]
                        << " dxB=" << polar_debug_.face_stencil_geom[side][3]
                        << " dxp=" << polar_debug_.face_stencil_geom[side][4]
                        << " dxm=" << polar_debug_.face_stencil_geom[side][5] << "\n";
              for (int nv=0; nv<2; ++nv) {
                std::cout << "plm " << side_name[side] << " " << var_name[nv]
                          << ": qm=" << polar_debug_.face_stencil_q[side][nv][0]
                          << " qc=" << polar_debug_.face_stencil_q[side][nv][1]
                          << " qp=" << polar_debug_.face_stencil_q[side][nv][2]
                          << " dwl=" << polar_debug_.face_stencil_calc[side][nv][0]
                          << " dwr=" << polar_debug_.face_stencil_calc[side][nv][1]
                          << " dqB=" << polar_debug_.face_stencil_calc[side][nv][2]
                          << " dqF=" << polar_debug_.face_stencil_calc[side][nv][3]
                          << " dwm=" << polar_debug_.face_stencil_calc[side][nv][4]
                          << " qminus=" << polar_debug_.face_stencil_calc[side][nv][5]
                          << " qplus=" << polar_debug_.face_stencil_calc[side][nv][6]
                          << "\n";
                std::cout << "plm_actual " << side_name[side] << " "
                          << var_name[nv]
                          << ": hydro_qplus="
                          << polar_debug_.face_plm_after_hydro[side][nv][0]
                          << " hydro_qminus="
                          << polar_debug_.face_plm_after_hydro[side][nv][1]
                          << " after_scalar_qplus="
                          << polar_debug_.face_plm_after_scalars[side][nv][0]
                          << " after_scalar_qminus="
                          << polar_debug_.face_plm_after_scalars[side][nv][1]
                          << " wc=" << polar_debug_.face_plm_scratch[side][nv][0]
                          << " dwl=" << polar_debug_.face_plm_scratch[side][nv][1]
                          << " dwr=" << polar_debug_.face_plm_scratch[side][nv][2]
                          << " dwm=" << polar_debug_.face_plm_scratch[side][nv][3]
                          << "\n";
              }
            }
          }

          std::cout << "POLAR_D_NEGATIVE_END" << std::endl;
#ifdef MPI_PARALLEL
          MPI_Abort(MPI_COMM_WORLD, 86);
#endif
          std::abort();
        }
      }

      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          u_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void Hydro::AddFluxDivergence_STS
//! \brief Adds diffusive flux divergence to weighted average of conservative variables
//! from previous step of time integration algorithm.  RKL2 registers are set to flux
//! divergence update if first stage of RKL2 STS.  Only a subset of ph->u indices are
//! integrated, dependent on which diffusive processes are integrated.
void Hydro::AddFluxDivergence_STS(const Real wght, int stage,
                                  AthenaArray<Real> &u_out,
                                  AthenaArray<Real> &fl_div_out,
                                  std::vector<int> idx_subset) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = flux[X1DIR];
  AthenaArray<Real> &x2flux = flux[X2DIR];
  AthenaArray<Real> &x3flux = flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k, j, is, ie+1, x1area);
      for (int n=0; n<NHYDRO; ++n) {
        if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) = (x1area(i+1) *x1flux(n,k,j,i+1) - x1area(i)*x1flux(n,k,j,i));
          }
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j  , is, ie, x2area   );
        pmb->pcoord->Face2Area(k, j+1, is, ie, x2area_p1);
        for (int n=0; n<NHYDRO; ++n) {
          if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              dflx(n,i) += (x2area_p1(i)*x2flux(n,k,j+1,i) - x2area(i)*x2flux(n,k,j,i));
            }
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  , j, is, ie, x3area   );
        pmb->pcoord->Face3Area(k+1, j, is, ie, x3area_p1);
        for (int n=0; n<NHYDRO; ++n) {
          if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
            for (int i=is; i<=ie; ++i) {
              dflx(n,i) += (x3area_p1(i)*x3flux(n,k+1,j,i) - x3area(i)*x3flux(n,k,j,i));
            }
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k, j, is, ie, vol);
      for (int n=0; n<NHYDRO; ++n) {
        if (std::binary_search(idx_subset.begin(), idx_subset.end(), n)) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            u_out(n,k,j,i) -= wght*dflx(n,i)/vol(i);
            if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
              fl_div_out(n,k,j,i) = -0.5*pmb->pmy_mesh->dt*dflx(n,i)/vol(i);
            }
          }
        }
      }
    }
  }
  return;
}

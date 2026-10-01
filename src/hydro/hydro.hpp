#ifndef HYDRO_HYDRO_HPP_
#define HYDRO_HYDRO_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro.hpp
//! \brief definitions for Hydro class

// C headers

// C++ headers
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/hydro/bvals_hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"
#include "srcterms/hydro_srcterms.hpp"

class MeshBlock;
class ParameterInput;

// TODO(felker): consider adding a struct FaceFlux w/ overloaded ctor in athena.hpp, or:
// using FaceFlux = AthenaArray<Real>[3];

//! \class Hydro
//! \brief hydro data and functions

class Hydro {
  friend class Field;
  friend class EquationOfState;
 public:
  Hydro(MeshBlock *pmb, ParameterInput *pin);

  // data
  // TODO(KGF): make this private, if possible
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Hydro

  // conserved and primitive variables
  AthenaArray<Real> u, w;        // time-integrator memory register #1
  AthenaArray<Real> u1, w1;      // time-integrator memory register #2
  AthenaArray<Real> u2;          // time-integrator memory register #3
  AthenaArray<Real> u0, fl_div; // rkl2 STS memory registers;
  // for the HL3D2 solver
  AthenaArray<Real> dvn, dvt;
  // (no more than MAX_NREGISTER allowed)

  AthenaArray<Real> flux[3];  // face-averaged flux vector

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_cons_, coarse_prim_;
  int refinement_idx{-1};

  // fourth-order intermediate quantities
  AthenaArray<Real> u_cc, w_cc;      // cell-centered approximations

  HydroBoundaryVariable hbvar;
  HydroSourceTerms hsrc;
  HydroDiffusion hdif;

  // functions
  void NewBlockTimeStep();    // computes new timestep on a MeshBlock
  void AddFluxDivergence(const Real wght, AthenaArray<Real> &u_out);
  void AddFluxDivergence_STS(const Real wght, int stage,
                             AthenaArray<Real> &u_out,
                             AthenaArray<Real> &fl_div_out,
                             std::vector<int> idx_subset);
  void CalculateFluxes(AthenaArray<Real> &w, FaceField &b,
                       AthenaArray<Real> &bcc, const int order);
  void CalculateFluxes_STS();
#if !MAGNETIC_FIELDS_ENABLED  // Hydro:
  void RiemannSolver(
      const int k, const int j, const int il, const int iu,
      const int ivx,
      AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
      const AthenaArray<Real> &dxw
#if EOS_SCALAR_INPUT_ENABLED
      , AthenaArray<Real> *rl = nullptr, AthenaArray<Real> *rr = nullptr
#endif
      );
#else  // MHD:
  void RiemannSolver(
      const int k, const int j, const int il, const int iu,
      const int ivx, const AthenaArray<Real> &bx,
      AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx,
      AthenaArray<Real> &ey, AthenaArray<Real> &ez,
      AthenaArray<Real> &wct, const AthenaArray<Real> &dxw
#if EOS_SCALAR_INPUT_ENABLED
      , AthenaArray<Real> *rl = nullptr, AthenaArray<Real> *rr = nullptr
#endif
      );
#endif
  void CalculateVelocityDifferences(const int k, const int j, const int il, const int iu,
    const int ivx, AthenaArray<Real> &dvn, AthenaArray<Real> &dvt);

 private:
  // Temporary, narrowly targeted diagnostics for the first polar cell of gid 0.
  // The CFL fields are saved when NewBlockTimeStep() runs so that a later failed
  // update can report the quantities that selected the timestep being used.
  struct PolarCflDebugData {
    bool cfl_valid;
    bool face_valid;
    int cfl_cycle;
    int face_cycle;
    Real cfl_time;
    Real face_time;
    Real r;
    Real theta;
    Real dtheta;
    Real center_width1;
    Real center_width2;
    Real g00;
    Real g02;
    Real g22;
    Real gi00;
    Real gi02;
    Real gi22;
    Real speed1;
    Real speed2;
    Real dt1_raw;
    Real dt2_raw;
    Real block_dt_raw;
    Real block_dt_scaled;
    Real rho;
    Real press;
    Real u1;
    Real u2;
    Real u3;
    Real volume;
    Real area1_l;
    Real area1_r;
    Real area2_l;
    Real area2_r;
    Real fv_dt1_raw;
    Real fv_dt2_raw;
    Real fv_dt_multid_scaled;
    int face_order;
    bool face_characteristic_projection;
    bool face_minmod;
    Real face_wl_pre[NWAVE];
    Real face_wr_pre[NWAVE];
    Real face_wl[NWAVE];
    Real face_wr[NWAVE];
    // Exact primitive PLM stencils at the outgoing x2 face.  The first index is
    // the cell on the left/right of the face, the second is rho/pressure, and
    // the third is the lower/center/upper x2 stencil point.
    Real face_stencil_q[2][2][3];
    // cf, cb, dxF, dxB, dxp, dxm for each of the two contributing cells.
    Real face_stencil_geom[2][6];
    // dwl, dwr, dqB, dqF, dwm, qminus, qplus for rho/pressure in each cell.
    Real face_stencil_calc[2][2][7];
    // Exact PLM outputs and reconstruction scratch values sampled immediately
    // after hydro reconstruction, plus the same outputs sampled again after the
    // EOS-scalar reconstruction call.  Indices are [j-js][rho/pressure]
    // [qplus/qminus] and [j-js][rho/pressure][wc/dwl/dwr/dwm].
    Real face_plm_after_hydro[2][2][2];
    Real face_plm_after_scalars[2][2][2];
    Real face_plm_scratch[2][2][4];
    Real lambda_p_l;
    Real lambda_m_l;
    Real lambda_p_r;
    Real lambda_m_r;
    Real hllc_wgas_l;
    Real hllc_wgas_r;
    Real hllc_cs2_l;
    Real hllc_cs2_r;
    Real hllc_cons_d_l;
    Real hllc_cons_d_r;
    Real hllc_lambda_l;
    Real hllc_lambda_r;
    Real hllc_lambda_star;
    Real hllc_pgas_star;
    Real hllc_contact_discriminant;
    Real hllc_cons_d_lstar;
    Real hllc_cons_d_rstar;
  } polar_debug_{};

  AthenaArray<Real> dt1_, dt2_, dt3_;  // scratch arrays used in NewTimeStep
  // scratch space used to compute fluxes
  AthenaArray<Real> dxw_;
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  // 2D
  AthenaArray<Real> wl_, wr_, wlb_;
#if EOS_SCALAR_INPUT_ENABLED
  AthenaArray<Real> rl_, rr_, rlb_;
#endif
  AthenaArray<Real> dflx_;
  AthenaArray<Real> bb_normal_;    // normal magnetic field, for (SR/GR)MHD
  AthenaArray<Real> lambdas_p_l_;  // most positive wavespeeds in left state
  AthenaArray<Real> lambdas_m_l_;  // most negative wavespeeds in left state
  AthenaArray<Real> lambdas_p_r_;  // most positive wavespeeds in right state
  AthenaArray<Real> lambdas_m_r_;  // most negative wavespeeds in right state
  // 2D GR
  AthenaArray<Real> g_, gi_;       // metric and inverse, for timesteps and some rsolvers
  AthenaArray<Real> cons_;         // conserved state, for some GR Riemann solvers

  // fourth-order hydro
  // 4D scratch arrays
  AthenaArray<Real> scr1_nkji_, scr2_nkji_;
  AthenaArray<Real> wl3d_, wr3d_;
  // 1D scratch arrays
  AthenaArray<Real> laplacian_l_fc_, laplacian_r_fc_;

  TimeStepFunc UserTimeStep_;

  void AddDiffusionFluxes();
  Real GetWeightForCT(Real dflx, Real rhol, Real rhor, Real dx, Real dt);
};
#endif // HYDRO_HYDRO_HPP_

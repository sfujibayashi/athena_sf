//========================================================================================
//! \file gr_flux_compare.cpp
//! \brief Direct face-flux test for GR Valencia Riemann solvers.
//========================================================================================

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../metric/metric.hpp"
#include "../parameter_input.hpp"

#if !GENERAL_RELATIVITY
#error "gr_flux_compare requires general relativity"
#endif
#if !VALENCIA_FORMULATION
#error "gr_flux_compare requires the Valencia formulation"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "gr_flux_compare is currently hydro-only"
#endif

namespace {

void SetSyntheticState(AthenaArray<Real> &w, int i,
                       Real rho, Real pgas,
                       Real uhat1, Real uhat2, Real uhat3,
                       const AthenaArray<Real> &g) {
  const Real g11 = g(I11, i);
  const Real g22 = g(I22, i);
  const Real g33 = g(I33, i);
  if (!(g11 > 0.0 && g22 > 0.0 && g33 > 0.0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_flux_compare\n"
        << "Non-positive spatial metric component: "
        << g11 << " " << g22 << " " << g33 << std::endl;
    ATHENA_ERROR(msg);
  }
  w(IDN, i) = rho;
  w(IPR, i) = pgas;
  w(IVX, i) = uhat1/std::sqrt(g11);
  w(IVY, i) = uhat2/std::sqrt(g22);
  w(IVZ, i) = uhat3/std::sqrt(g33);
}

void RunFluxTest(MeshBlock *pmb, ParameterInput *pin, int dir) {
  if (dir == 2 && pmb->block_size.nx2 <= 1) return;
  if (dir == 3 && pmb->block_size.nx3 <= 1) return;

  int i = pmb->is + (pmb->ie-pmb->is)/2;
  int j = pmb->js + (pmb->je-pmb->js)/2;
  int k = pmb->ks + (pmb->ke-pmb->ks)/2;
  int ivx = IVX;

  if (dir == 1) {
    i = pmb->is + (pmb->ie-pmb->is+1)/2;
    ivx = IVX;
  } else if (dir == 2) {
    j = pmb->js + (pmb->je-pmb->js+1)/2;
    ivx = IVY;
  } else if (dir == 3) {
    k = pmb->ks + (pmb->ke-pmb->ks+1)/2;
    ivx = IVZ;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_flux_compare: dir must be 1,2,3" << std::endl;
    ATHENA_ERROR(msg);
  }

  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, pmb->ncells1+1);
  gi.NewAthenaArray(NMETRIC, pmb->ncells1+1);
  g.ZeroClear();
  gi.ZeroClear();

  if (dir == 1) pmb->pmetric->Face1Metric(k,j,i,i,g,gi);
  if (dir == 2) pmb->pmetric->Face2Metric(k,j,i,i,g,gi);
  if (dir == 3) pmb->pmetric->Face3Metric(k,j,i,i,g,gi);

  const Real rho_l = pin->GetOrAddReal("problem", "rho_l", 1.0);
  const Real p_l   = pin->GetOrAddReal("problem", "p_l",   0.3);
  const Real uh1_l = pin->GetOrAddReal("problem", "uh1_l", 0.4);
  const Real uh2_l = pin->GetOrAddReal("problem", "uh2_l", 0.2);
  const Real uh3_l = pin->GetOrAddReal("problem", "uh3_l",-0.1);
  const Real rho_r = pin->GetOrAddReal("problem", "rho_r", 0.7);
  const Real p_r   = pin->GetOrAddReal("problem", "p_r",   0.15);
  const Real uh1_r = pin->GetOrAddReal("problem", "uh1_r",-0.2);
  const Real uh2_r = pin->GetOrAddReal("problem", "uh2_r", 0.1);
  const Real uh3_r = pin->GetOrAddReal("problem", "uh3_r", 0.25);

  AthenaArray<Real> wl, wr, flx, dxw;
  wl.NewAthenaArray(NWAVE, pmb->ncells1+1);
  wr.NewAthenaArray(NWAVE, pmb->ncells1+1);
  flx.NewAthenaArray(NHYDRO, pmb->ncells3+1, pmb->ncells2+1, pmb->ncells1+1);
  dxw.NewAthenaArray(pmb->ncells1+1);
  wl.ZeroClear(); wr.ZeroClear(); flx.ZeroClear(); dxw.ZeroClear();

  SetSyntheticState(wl,i,rho_l,p_l,uh1_l,uh2_l,uh3_l,g);
  SetSyntheticState(wr,i,rho_r,p_r,uh1_r,uh2_r,uh3_r,g);

  // Call the actual compile-time-selected Riemann solver.
  pmb->phydro->RiemannSolver(k,j,i,i,ivx,wl,wr,flx,dxw);

  const Real alpha = std::sqrt(-1.0/gi(I00,i));
  std::cout << std::scientific << std::setprecision(17)
            << "FLUX_TEST"
            << " dir=" << dir
            << " i=" << i << " j=" << j << " k=" << k
            << " alpha=" << alpha
            << " g00=" << g(I00,i)
            << " g11=" << g(I11,i)
            << " g22=" << g(I22,i)
            << " g33=" << g(I33,i)
            << " FD=" << flx(IDN,k,j,i)
            << " Ftau=" << flx(IEN,k,j,i)
            << " FS1=" << flx(IM1,k,j,i)
            << " FS2=" << flx(IM2,k,j,i)
            << " FS3=" << flx(IM3,k,j,i)
            << std::endl;
}

} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  AllocateRealUserMeshDataField(2);
  ruser_mesh_data[0].NewAthenaArray(1);
  ruser_mesh_data[1].NewAthenaArray(1);
  ruser_mesh_data[0](0) = pin->GetOrAddReal("coord", "m", 1.0);
  ruser_mesh_data[1](0) = pin->GetOrAddReal("coord", "j", 0.0);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(ncells1+1); // Phi_face1
  ruser_meshblock_data[1].NewAthenaArray(ncells1+1); // Mgrav_face1
  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Put the actual grid at the EOS atmosphere floor.  The RGPS atmosphere
  // subtraction then leaves essentially only the specified BH metric.
  for (int k=0; k<ncells3; ++k) {
    for (int j=0; j<ncells2; ++j) {
      for (int i=0; i<ncells1; ++i) {
        const Real r = pcoord->x1v(i);
        const Real rho = peos->GetDensityFloor(r);
        const Real pgas = peos->GetPressureFloor(r);
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = 0.0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = 0.0;
      }
    }
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  if (gid != 0) return;
  static bool done = false;
  if (done) return;
  done = true;

  RunFluxTest(this,pin,1);
  RunFluxTest(this,pin,2);
  RunFluxTest(this,pin,3);
}

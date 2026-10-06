//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_round_trip.cpp
//! \brief Primitive -> conserved -> primitive round-trip test for GR hydro.
//!
//! The input velocity
//! components uhat1, uhat2, uhat3 are orthonormal-frame components of the projected spatial
//! four-velocity \tilde{u}^{\hat{i}} = W v^{\hat{i}}.  They are converted to coordinate
//! components using the (diagonal) monopole spatial metric before calling the EOS routines.
//========================================================================================

// C++ headers
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !GENERAL_RELATIVITY
#error "eos_round_trip.cpp requires GENERAL_RELATIVITY"
#endif

#if MAGNETIC_FIELDS_ENABLED
#error "eos_round_trip.cpp is intended for pure hydrodynamics"
#endif

namespace {

inline Real ScaledError(Real x, Real xref) {
  return std::abs(x - xref) / std::max(static_cast<Real>(1.0), std::abs(xref));
}

}  // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  //
  AllocateRealUserMeshDataField(2);

  // BH mass
  ruser_mesh_data[0].NewAthenaArray(1);
  // BH spin
  ruser_mesh_data[1].NewAthenaArray(1);

  ruser_mesh_data[0](0) = pin->GetOrAddReal("coord", "m", 0.0);
  ruser_mesh_data[1](0) = pin->GetOrAddReal("coord", "j", 0.0);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  AllocateRealUserMeshBlockDataField(2);
  
  // [0] for Psi_face1
  ruser_meshblock_data[0].NewAthenaArray(ncells1 + 1);
  // [1] for delta_m_face1
  ruser_meshblock_data[1].NewAthenaArray(ncells1 + 1);

  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();

  
}
//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Test primitive -> conserved -> primitive conversion.

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Test state.  The velocity variables are orthonormal-frame components of
  // \tilde{u}^{\hat{i}} = W v^{\hat{i}}, not ordinary three-velocities.
  const Real rho0   = pin->GetOrAddReal("problem", "rho",   1.0);
  const Real pgas0  = pin->GetOrAddReal("problem", "pgas", 0.1);
  const Real uhat1  = pin->GetOrAddReal("problem", "uhat1", 0.2);
  const Real uhat2  = pin->GetOrAddReal("problem", "uhat2", 0.1);
  const Real uhat3  = pin->GetOrAddReal("problem", "uhat3", 0.05);
  const Real tol    = pin->GetOrAddReal("problem", "round_trip_tol", 1.0e-12);
  const bool abort_on_fail =
      pin->GetOrAddBoolean("problem", "abort_on_fail", true);

  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);

  // -----------------------------------------------------------------------------
  // 1. Set input primitives and save an untouched copy in w1.
  //
  // For the monopole metric used here the spatial metric is diagonal, so
  //   uhat1 = sqrt(g11) * uu1,
  //   uhat2 = sqrt(g22) * uu2,
  //   uhat3 = sqrt(g33) * uu3.
  // -----------------------------------------------------------------------------
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);

      for (int i=is; i<=ie; ++i) {
        const Real h1 = std::sqrt(g(I11,i));
        const Real h2 = std::sqrt(g(I22,i));
        const Real h3 = std::sqrt(g(I33,i));

        phydro->w(IDN,k,j,i) = rho0;
        phydro->w(IPR,k,j,i) = pgas0;
        phydro->w(IVX,k,j,i) = uhat1/h1;
        phydro->w(IVY,k,j,i) = uhat2/h2;
        phydro->w(IVZ,k,j,i) = uhat3/h3;

        for (int n=0; n<NHYDRO; ++n) {
          phydro->w1(n,k,j,i) = phydro->w(n,k,j,i);
        }
      }
    }
  }

  // -----------------------------------------------------------------------------
  // 2. Primitive -> conserved.
  // -----------------------------------------------------------------------------
  AthenaArray<Real> dummy_bcc;
  dummy_bcc.NewAthenaArray(3, ncells3, ncells2, ncells1);
  dummy_bcc.ZeroClear();
  FaceField dummy_b;

  peos->PrimitiveToConserved(phydro->w, dummy_bcc, phydro->u, pcoord,
                             is, ie, js, je, ks, ke);

  // -----------------------------------------------------------------------------
  // 3. Conserved -> primitive.
  //    w1 is used both as the pressure guess and as the reference state.
  // -----------------------------------------------------------------------------
  peos->ConservedToPrimitive(phydro->u, phydro->w1, dummy_b,
                             phydro->w, dummy_bcc, pcoord,
                             is, ie, js, je, ks, ke);

  // -----------------------------------------------------------------------------
  // 4. Compare recovered primitives with the input primitives.
  // -----------------------------------------------------------------------------
  Real max_err = 0.0;
  Real max_err_var[NHYDRO];
  for (int n=0; n<NHYDRO; ++n) max_err_var[n] = 0.0;

  int imax = is, jmax = js, kmax = ks, nmax = 0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<NHYDRO; ++n) {
          const Real err = ScaledError(phydro->w(n,k,j,i), phydro->w1(n,k,j,i));
          max_err_var[n] = std::max(max_err_var[n], err);
          if (err > max_err) {
            max_err = err;
            imax = i;
            jmax = j;
            kmax = k;
            nmax = n;
          }
        }
      }
    }
  }

  std::cout << std::setprecision(17)
            << "  rho   : " << rho0 << "\n"
            << "  p     : " << pgas0 << "\n"
            << "  uhat1 : " << uhat1 << "\n"
            << "  uhat2 : " << uhat2 << "\n"
            << "  uhat3 : " << uhat3 << std::endl;


  std::cout << std::setprecision(17)
            << "P2C->C2P round-trip: gid=" << gid
            << " max_err=" << max_err
            << " at (n,k,j,i)=(" << nmax << "," << kmax << ","
            << jmax << "," << imax << ")\n"
            << "  rho : " << max_err_var[IDN] << "\n"
            << "  p   : " << max_err_var[IPR] << "\n"
            << "  uu1 : " << max_err_var[IVX] << "\n"
            << "  uu2 : " << max_err_var[IVY] << "\n"
            << "  uu3 : " << max_err_var[IVZ] << std::endl;

  if (abort_on_fail && max_err > tol) {
    std::stringstream msg;
    msg << "### FATAL ERROR in valencia_round_trip.cpp\n"
        << "P2C->C2P round-trip failed on gid=" << gid << "\n"
        << "max_err=" << max_err << " > tolerance=" << tol << "\n";
    ATHENA_ERROR(msg);
  }
  std::abort();
  return;
}

void MeshBlock::InitializeAtmosphere(ParameterInput *pin) {
  return;
}

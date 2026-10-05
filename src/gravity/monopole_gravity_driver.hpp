#ifndef GRAVITY_MONOPOLE_GRAVITY_DRIVER_HPP_
#define GRAVITY_MONOPOLE_GRAVITY_DRIVER_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monopole_gravity.hpp
//! \brief defines MonopoleGravityDriver class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "dynamic_metric_driver.hpp"

class Mesh;
class ParameterInput;

//! \class MonopoleGravityDriver
//! \brief Constructs global radial monopole gravity profiles
class MonopoleGravityDriver : public DynamicMetricDriver {
public:
  MonopoleGravityDriver(Mesh *pm, ParameterInput *pin);
  ~MonopoleGravityDriver() override;
  
  void Update() override;

  Real bh_mass_prev_, bh_spin_prev_;
  Real bh_mass_pending_, bh_spin_pending_;
  Real mdot_bh_, angdot_bh_;

  Real BlackHoleMassAccretionRate() const;
  void UpdateBlackHoleMass(int stage);

  Real GetBlackHoleMass() const;
  Real GetBlackHoleSpin() const;
  Real GetBlackHoleMassAccretionRate() const override;

  void CellMetricRadialDerivatives(
         MeshBlock *pmb,
         const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &d1_g00,
	 AthenaArray<Real> &d1_g11) const override;


  AthenaArray<Real>& PsiFace1(MeshBlock *pmb);
  const AthenaArray<Real>& PsiFace1(MeshBlock *pmb) const;
  
  AthenaArray<Real>& DeltaMFace1(MeshBlock *pmb);
  const AthenaArray<Real>& DeltaMFace1(MeshBlock *pmb) const;

  Real CellPsi(MeshBlock *pmb, int i) const;
  Real CellDeltaM(MeshBlock *pmb, int i) const;

private:
  Mesh *pmy_mesh_;

  int nr_;

  AthenaArray<Real> dm_shell_global_;
  AthenaArray<Real> delta_m_face_global_;
  AthenaArray<Real> Psi_face_global_;

  Real& BlackHoleMassStorage();
  const Real& BlackHoleMassStorage() const;

  Real& BlackHoleSpinStorage();
  const Real& BlackHoleSpinStorage() const;

  MonopoleGravity gravity_model_;
};

#endif // GRAVITY_MONOPOLE_GRAVITY_DRIVER_HPP_

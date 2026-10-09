#ifndef GRAVITY_RGPS_GRAVITY_DRIVER_HPP_
#define GRAVITY_RGPS_GRAVITY_DRIVER_HPP_
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
#include "rgps_gravity.hpp"

class Mesh;
class ParameterInput;

//! \class RGPSGravityDriver
//! \brief Constructs global RGPS gravity profiles
class RGPSGravityDriver : public DynamicMetricDriver {
public:
  RGPSGravityDriver(Mesh *pm, ParameterInput *pin);
  ~RGPSGravityDriver() override;
  
  void InitializeRadialGrid();
  void Update() override;
  void UpdateAfterConservedToPrimitive() override;

  Real bh_mass_prev_, bh_spin_prev_;
  Real bh_mass_pending_, bh_spin_pending_;
  Real mdot_bh_, angdot_bh_;

  Real BlackHoleMassAccretionRate() const;
  void UpdateBlackHoleMass(int stage);

  Real GetBlackHoleMass() const override;
  Real GetBlackHoleSpin() const;
  Real GetBlackHoleMassAccretionRate() const override;

  void CellMetricRadialDerivatives(
         MeshBlock *pmb,
         const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &d1_g00,
	 AthenaArray<Real> &d1_g11) const override;

  void Face1Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const override;

  void Face2Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const override;

  void Face3Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const override;

  void CellMetric(
       MeshBlock *pmb,
       const int k, const int j,
       const int il, const int iu,
       AthenaArray<Real> &g,
       AthenaArray<Real> &g_inv) const override;

  void ConstructCellCovariantMetric(
         MeshBlock *pmb,
	 int k, int j, int i,
	 Real &g00, Real &g01, Real &g02, Real &g03,
	 Real &g11, Real &g12, Real &g13,
	 Real &g22, Real &g23, Real &g33) const override;
  
  Real SqrtMinusG(
         MeshBlock *pmb,
	 const int k, const int j, const int i) const override;

  Real CellDensitizationFactor(
        MeshBlock *pmb, int k, int j, int i) const override;

  Real Face1DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const override;
  Real Face2DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const override;
  Real Face3DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const override;
  

  Real GetEnclosedMassAtInnerBoundary(
        MeshBlock *pmb) const override;


  AthenaArray<Real>& PhiFace1(MeshBlock *pmb);
  const AthenaArray<Real>& PhiFace1(MeshBlock *pmb) const;
  
  AthenaArray<Real>& MgravFace1(MeshBlock *pmb);
  const AthenaArray<Real>& MgravFace1(MeshBlock *pmb) const ;

  Real CellPhi(MeshBlock *pmb, int i) const;
  Real CellMgrav(MeshBlock *pmb, int i) const;

  // output
  int NumModelOutputVariables() const override;

  const char *ModelOutputVariableName(int n) const override;

  Real ModelOutputVariable(
    MeshBlock *pmb,
    int n,
    int k, int j, int i) const override;


private:
  Mesh *pmy_mesh_;

  int nr_;

  AthenaArray<Real> calE_shell_global_;
  AthenaArray<Real> Srr_shell_global_;
  AthenaArray<Real> mgrav_face_global_;
  AthenaArray<Real> Phi_face_global_;

  AthenaArray<Real> r_cell_global_;
  AthenaArray<Real> r_face_global_;

  AthenaArray<Real> vol_;


  Real& BlackHoleMassStorage();
  const Real& BlackHoleMassStorage() const;

  Real& BlackHoleSpinStorage();
  const Real& BlackHoleSpinStorage() const;

  RGPSGravity gravity_model_;
};

#endif // GRAVITY_RGPS_GRAVITY_DRIVER_HPP_

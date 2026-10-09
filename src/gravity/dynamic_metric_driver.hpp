#ifndef GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_
#define GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

#include "../athena.hpp" // Real
#include "../athena_arrays.hpp"

class Mesh;
class MeshBlock;
class ParameterInput;


class DynamicMetricDriver {
public:
  virtual ~DynamicMetricDriver() = default;
  
  virtual void InitializeFromPrimitive() = 0;
  virtual void InitializeFromRestart() = 0;
  virtual void UpdateBeforeCons2Prim(int stage) = 0;
  virtual void UpdateAfterCons2Prim() = 0;
  virtual void UpdateBlackHoleMass(int stage) = 0;

  virtual Real GetBlackHoleMassAccretionRate() const = 0;

  virtual void CellMetricRadialDerivatives(
         MeshBlock *pmb,
         const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &d1_g00,
	 AthenaArray<Real> &d1_g11) const = 0;

  virtual void Face1Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const = 0;

  virtual void Face2Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const = 0;
  
  virtual void Face3Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const = 0;

  virtual void CellMetric(
         MeshBlock *pmb,
	 const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &g,
	 AthenaArray<Real> &g_inv) const = 0;

  virtual void ConstructCellCovariantMetric(
         MeshBlock *pmb,
	 int k, int j, int i,
	 Real &g00, Real &g01, Real &g02, Real &g03,
	 Real &g11, Real &g12, Real &g13,
	 Real &g22, Real &g23, Real &g33) const = 0;

  virtual Real SqrtMinusG(
         MeshBlock *pmb,
	 const int k, const int j, const int i) const = 0;
  
  virtual Real CellDensitizationFactor(
        MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face1DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face2DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face3DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;

  virtual Real CellSpatialDensitizationFactor(
        MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face1SpatialDensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face2SpatialDensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;
  virtual Real Face3SpatialDensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const = 0;

  virtual void CellExtrinsicCurvature(
         MeshBlock *pmb,
	 const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &k11, AthenaArray<Real> &k22, AthenaArray<Real> &k33) const = 0;

  virtual Real GetEnclosedMassAtInnerBoundary(
        MeshBlock *pmb) const = 0;

  virtual Real GetBlackHoleMass() const = 0;


  // output
  virtual int NumModelOutputVariables() const = 0;
  
  virtual const char *ModelOutputVariableName(int n) const = 0;
  
  virtual Real ModelOutputVariable(
    MeshBlock *pmb,
    int n,
    int k, int j, int i) const = 0;

};

DynamicMetricDriver *CreateDynamicMetricDriver(Mesh *pm, ParameterInput *pin);

#endif // GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

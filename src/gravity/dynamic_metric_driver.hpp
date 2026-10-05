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
  
  virtual void Update() = 0;
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
  
  virtual void ConstructCellCovariantMetric(
         MeshBlock *pmb,
	 int k, int j, int i,
	 Real &g00, Real &g01, Real &g02, Real &g03,
	 Real &g11, Real &g12, Real &g13,
	 Real &g22, Real &g23, Real &g33) const = 0;

  
  
};

DynamicMetricDriver *CreateDynamicMetricDriver(Mesh *pm, ParameterInput *pin);

#endif // GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

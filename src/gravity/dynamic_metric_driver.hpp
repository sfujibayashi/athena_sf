#ifndef GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_
#define GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

#include "../athena.hpp" // Real
#include "../athena_arrays.hpp"

class Mesh;
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
  
};

DynamicMetricDriver *CreateDynamicMetricDriver(Mesh *pm, ParameterInput *pin);

#endif // GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

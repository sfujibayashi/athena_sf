#ifndef GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_
#define GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

#include "../athena.hpp" // Real

class DynamicMetricDriver {
 public:
  virtual ~DynamicMetricDriver() = default;

  virtual void Update() = 0;
  virtual void UpdateBlackHoleMass(int stage) = 0;

  virtual Real GetBlackHoleMassAccretionRate() const = 0;
  
};


#endif // GRAVITY_DYNAMIC_METRIC_DRIVER_HPP_

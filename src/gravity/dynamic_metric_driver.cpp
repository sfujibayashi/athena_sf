
// C++ headers
#include <sstream>

// Athena++ headers

#include "dynamic_metric_driver.hpp"
#include "monopole_gravity.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

DynamicMetricDriver *CreateDynamicMetricDriver(Mesh *pm, ParameterInput *pin){
  return new MonopoleGravityDriver(pm, pin);
}

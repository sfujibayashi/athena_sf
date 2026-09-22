#include "metric.hpp"

#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

Metric::Metric(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block(pmb) {
  const int nc1 = pmb->ncells1;
  const int nc2 = pmb->ncells2;
  const int nc3 = pmb->ncells3;

  alpha.NewAthenaArray(nc3, nc2, nc1);
  beta.NewAthenaArray(3, nc3, nc2, nc1);
  gamma.NewAthenaArray(N_GAMMA, nc3, nc2, nc1);
}

Metric::~Metric() {
}

void Metric::Update(Real time) {
  // For now empty.
  // Later:
  // 1. construct gravity source
  // 2. solve self-gravity equation
  // 3. construct alpha, beta^i, gamma_ij
}

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

  // initialize with Minkowski metric
  for(int k=0; k<nc3; ++k){
    for(int j=0; j<nc2; ++j){
      for(int i=0; i<nc1; ++i){
        alpha(k,j,i) = 1.0;
        beta(0,k,j,i) = 0.0;
        beta(1,k,j,i) = 0.0;
        beta(2,k,j,i) = 0.0;

        gamma(I_G11,k,j,i) = 1.0;
        gamma(I_G22,k,j,i) = 1.0;
        gamma(I_G33,k,j,i) = 1.0;

        gamma(I_G12,k,j,i) = 0.0;
        gamma(I_G13,k,j,i) = 0.0;
        gamma(I_G23,k,j,i) = 0.0;
      }
    }
  }
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

#include <cmath>

#include "metric.hpp"

#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

Metric::Metric(MeshBlock *pmb, ParameterInput *pin)
    : pmy_block(pmb),
      bh_mass_(pmb->pcoord->GetMass())  {
  const int nc1 = pmb->ncells1;
  const int nc2 = pmb->ncells2;
  const int nc3 = pmb->ncells3;

  alpha.NewAthenaArray(nc3, nc2, nc1);
  beta.NewAthenaArray(3, nc3, nc2, nc1);
  gamma.NewAthenaArray(N_GAMMA, nc3, nc2, nc1);

  // for now
  Update(0.0);
}

Metric::~Metric() {
}

void Metric::Update(Real time) {
  Coordinates *pcoord = pmy_block->pcoord;
  
  const int nc1 = pmy_block->ncells1;
  const int nc2 = pmy_block->ncells2;
  const int nc3 = pmy_block->ncells3;
  
  for (int k=0; k<nc3; ++k) {
    for (int j=0; j<nc2; ++j) {
      const Real theta = pcoord->x2v(j);
      const Real sin_theta = std::sin(theta);
      
      for (int i=0; i<nc1; ++i) {
        const Real r = pcoord->x1v(i);
        const Real alpha_sq = 1.0 - 2.0*bh_mass_/r;
        
        alpha(k,j,i) = std::sqrt(alpha_sq);
        
        beta(0,k,j,i) = 0.0;
        beta(1,k,j,i) = 0.0;
        beta(2,k,j,i) = 0.0;
        
        gamma(I_G11,k,j,i) = 1.0/alpha_sq;
        gamma(I_G12,k,j,i) = 0.0;
        gamma(I_G13,k,j,i) = 0.0;
        gamma(I_G22,k,j,i) = r*r;
        gamma(I_G23,k,j,i) = 0.0;
        gamma(I_G33,k,j,i) = r*r*sin_theta*sin_theta;
      }
    }
  }
  
  // 1. construct gravity source
  // 2. solve self-gravity equation
  // 3. construct alpha, beta^i, gamma_ij
}

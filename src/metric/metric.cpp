// C++ headers
#include <cmath>
#include <sstream>

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

Real Metric::DetSpatialMetric(Real g11, Real g12, Real g13,
                              Real g22, Real g23, Real g33) const {
  const Real det =
    g11 * (g22*g33 - g23*g23)
    - g12 * (g12*g33 - g13*g23)
    + g13 * (g12*g23 - g13*g22);
    
  if (det <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in InvertSpatialMetric\n"
        << "Non-positive spatial metric determinant: det(gamma) = "
        << det << std::endl;
    ATHENA_ERROR(msg);
  }
  return det;
}


void Metric::InvertSpatialMetric(Real g11, Real g12, Real g13,
                                 Real g22, Real g23, Real g33,
                                 Real &gi11, Real &gi12, Real &gi13,
                                 Real &gi22, Real &gi23, Real &gi33) const {
  
  const Real det = DetSpatialMetric(g11,g12,g13,g22,g23,g33);
  
  const Real inv_det = 1.0 / det;
    
  gi11 =  (g22*g33 - g23*g23) * inv_det;
  gi12 =  (g13*g23 - g12*g33) * inv_det;
  gi13 =  (g12*g23 - g13*g22) * inv_det;
    
  gi22 =  (g11*g33 - g13*g13) * inv_det;
  gi23 =  (g12*g13 - g11*g23) * inv_det;
    
  gi33 =  (g11*g22 - g12*g12) * inv_det;
}
  
void Metric::CellMetric(const int k, const int j, const int il, const int iu,
                AthenaArray<Real> &g, AthenaArray<Real> &g_inv){
  for(int i=il; i<=iu; ++i){
    const Real a = alpha(k,j,i);

    // beta^i
    const Real b1 = beta(0,k,j,i);
    const Real b2 = beta(1,k,j,i);
    const Real b3 = beta(2,k,j,i);
      
    // g_ij = gamma_ij
    const Real g11 = gamma(I_G11,k,j,i);
    const Real g12 = gamma(I_G12,k,j,i);
    const Real g13 = gamma(I_G13,k,j,i);
    const Real g22 = gamma(I_G22,k,j,i);
    const Real g23 = gamma(I_G23,k,j,i);
    const Real g33 = gamma(I_G33,k,j,i);

    // g_0i = gamma_ij beta^j
    const Real g01 = g11*b1 + g12*b2 + g13*b3;
    const Real g02 = g12*b1 + g22*b2 + g23*b3;
    const Real g03 = g13*b1 + g23*b2 + g33*b3;
      
    // gamma_ij beta^j beta^i = g_0i beta^i
    const Real g00 = -a*a + g01*b1 + g02*b2 + g03*b3;

    g(I00,i) = g00;
    g(I01,i) = g01;
    g(I02,i) = g02;
    g(I03,i) = g03;
    g(I11,i) = g11;
    g(I12,i) = g12;
    g(I13,i) = g13;
    g(I22,i) = g22;
    g(I23,i) = g23;
    g(I33,i) = g33;
      
    Real gi11, gi12, gi13, gi22, gi23, gi33;
    InvertSpatialMetric(g11, g12, g13, g22, g23, g33,
                        gi11, gi12, gi13, gi22, gi23, gi33);
    const Real a2i = 1.0/(a*a);
    g_inv(I00,i) = -a2i;
    g_inv(I01,i) = b1*a2i;
    g_inv(I02,i) = b2*a2i;
    g_inv(I03,i) = b3*a2i;
    g_inv(I11,i) = gi11 - b1*b1*a2i;
    g_inv(I12,i) = gi12 - b1*b2*a2i;
    g_inv(I13,i) = gi13 - b1*b3*a2i;
    g_inv(I22,i) = gi22 - b2*b2*a2i;
    g_inv(I23,i) = gi23 - b2*b3*a2i;
    g_inv(I33,i) = gi33 - b3*b3*a2i;
      
  }
}

Real Metric::SqrtMinusG(int k, int j, int i) const {
  
  // lapse
  const Real a = alpha(k,j,i);
  
  // g_ij = gamma_ij
  const Real g11 = gamma(I_G11,k,j,i);
  const Real g12 = gamma(I_G12,k,j,i);
  const Real g13 = gamma(I_G13,k,j,i);
  const Real g22 = gamma(I_G22,k,j,i);
  const Real g23 = gamma(I_G23,k,j,i);
  const Real g33 = gamma(I_G33,k,j,i);
  
  const Real detgamma = DetSpatialMetric(g11,g12,g13,g22,g23,g33);
  
  return a*std::sqrt(detgamma);
}

Real Metric::DensitizationFactor(int k, int j, int i) const {
  const Real sqrt_minus_g = SqrtMinusG(k,j,i);
  const Real r = pmy_block->pcoord->x1v(i);
  const Real theta = pmy_block->pcoord->x2v(j);
  return sqrt_minus_g/(r*r*std::sin(theta));
}

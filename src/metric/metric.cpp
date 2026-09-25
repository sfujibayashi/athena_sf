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

  Psi_.NewAthenaArray(nc1);
  delta_m_.NewAthenaArray(nc1);

  Psi_.ZeroClear();
  delta_m_.ZeroClear();
}

Metric::~Metric() {
}

// delta_m_ and Psi_ are derived from fluid distribution.
void Metric::Update(Real time) {
  Coordinates *pcoord = pmy_block->pcoord;
  
  // const int nc1 = pmy_block->ncells1;
  // const int nc2 = pmy_block->ncells2;
  // const int nc3 = pmy_block->ncells3;
  
  // for (int k=0; k<nc3; ++k) {
  //   for (int j=0; j<nc2; ++j) {
  //     for (int i=0; i<nc1; ++i) {
  //     }
  //   }
  // }
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
  
  Coordinates *pcoord = pmy_block->pcoord;
  
  const Real theta = pcoord->x2v(j);
  for(int i=il; i<=iu; ++i){
    const Real r = pcoord->x1v(i);
    
    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;
    ConstructCellCovariantMetric(k,j,i,
      g00, g01, g02, g03,
      g11, g12, g13, g22, g23, g33);
    
    g(I00) = g00;
    g(I01) = g01;
    g(I02) = g02;
    g(I03) = g03;
    g(I11) = g11;
    g(I12) = g12;
    g(I13) = g13;
    g(I22) = g22;
    g(I23) = g23;
    g(I33) = g33;
    
    InvertMetric(i, g, g_inv);
  }
}

void Metric::Face1Metric(const int k, const int j, const int il, const int iu,
                                AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  Coordinates *pcoord = pmy_block->pcoord;
  
  const Real theta = pcoord->x2v(j);
  
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    const Real r = pcoord->x1f(i);
    ConstructBackgroundMetric(r, theta, i, g);
    
    const Real Psi = 0.0;
    const Real delta_m = 0.0;
    AddSelfGravityPerturbation(r, Psi, delta_m, i, g);
    
    InvertMetric(i, g, g_inv);
  }
  return;
}

void Metric::Face2Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  Coordinates *pcoord = pmy_block->pcoord;
  
  const Real theta = pcoord->x2f(j);
  
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

    const Real r = pcoord->x1v(i);

    ConstructBackgroundMetric(r, theta, i, g);
    
    const Real Psi = 0.0;
    const Real delta_m = 0.0;
    AddSelfGravityPerturbation(r, Psi, delta_m, i, g);
    
    InvertMetric(i, g, g_inv);
  }
  return;
}

void Metric::Face3Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  Coordinates *pcoord = pmy_block->pcoord;
  
  const Real theta = pcoord->x2v(j);

  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

    const Real r = pcoord->x1v(i);

    ConstructBackgroundMetric(r, theta, i, g);
    
    const Real Psi = 0.0;
    const Real delta_m = 0.0;
    AddSelfGravityPerturbation(r, Psi, delta_m, i, g);
    
    InvertMetric(i, g, g_inv);

  }
  return;
}

Real Metric::SqrtMinusG(int k, int j, int i) const {
  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  ConstructCellCovariantMetric(k, j, i,
      g00, g01, g02, g03,
      g11, g12, g13, g22, g23, g33);
  
  const Real detgamma = DetSpatialMetric(g11, g12, g13, g22, g23, g33);

  // beta_i = g_0i
  Real gi11, gi12, gi13, gi22, gi23, gi33;
  InvertSpatialMetric(g11, g12, g13, g22, g23, g33,
                      gi11, gi12, gi13, gi22, gi23, gi33);

  const Real beta1 = gi11*g01 + gi12*g02 + gi13*g03;
  const Real beta2 = gi12*g01 + gi22*g02 + gi23*g03;
  const Real beta3 = gi13*g01 + gi23*g02 + gi33*g03;


  const Real alpha_sq = -g00 + g01*beta1 + g02*beta2 + g03*beta3;

  return std::sqrt(alpha_sq*detgamma);
}

Real Metric::DensitizationFactor(int k, int j, int i) const {
  const Real sqrt_minus_g = SqrtMinusG(k,j,i);
  const Real r = pmy_block->pcoord->x1v(i);
  const Real theta = pmy_block->pcoord->x2v(j);
  return sqrt_minus_g/(r*r*std::sin(theta));
}


void Metric::SetBlackHoleMass(Real mass){
  bh_mass_ = mass;
}


Real Metric::GetBlackHoleMass() const {
  return bh_mass_;
}



// constructor of background metric g_munu
void Metric::ConstructBackgroundMetric(
    Real r, Real theta,
    int i,
    AthenaArray<Real> &g) const{
  
  const Real r_sq = SQR(r);
  const Real f = 1.0 - 2.0*bh_mass_/r;
  const Real sintheta = std::sin(theta);
  const Real sin2theta = sintheta*sintheta;

  g(I00,i) = -f;
  g(I01,i) = 0.0;
  g(I02,i) = 0.0;
  g(I03,i) = 0.0;
  g(I11,i) = 1.0/f;
  g(I12,i) = 0.0;
  g(I13,i) = 0.0;
  g(I22,i) = r_sq;
  g(I23,i) = 0.0;
  g(I33,i) = r_sq*sin2theta;

}

void Metric::AddSelfGravityPerturbation(
    Real r, Real Psi, Real delta_m,
    int i,
    AthenaArray<Real> &g) const {

  const Real f = 1.0 - 2.0*bh_mass_/r;

  const Real h00 = 2.0*delta_m/r + 2.0*f*Psi;
  const Real h11 = 2.0*delta_m / (r*f*f);
  g(I00,i) += h00;
  g(I11,i) += h11;

}


void Metric::InvertMetric(
    int i,
    const AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv) const {

  // spatial metric
  const Real g11 = g(I11,i);
  const Real g12 = g(I12,i);
  const Real g13 = g(I13,i);
  const Real g22 = g(I22,i);
  const Real g23 = g(I23,i);
  const Real g33 = g(I33,i);

  Real gi11, gi12, gi13, gi22, gi23, gi33;

  InvertSpatialMetric(
      g11, g12, g13,
      g22, g23, g33,
      gi11, gi12, gi13,
      gi22, gi23, gi33);

  // beta_i = g_0i
  const Real b_1 = g(I01,i);
  const Real b_2 = g(I02,i);
  const Real b_3 = g(I03,i);

  // beta^i = gamma^{ij} beta_j
  const Real b1 = gi11*b_1 + gi12*b_2 + gi13*b_3;
  const Real b2 = gi12*b_1 + gi22*b_2 + gi23*b_3;
  const Real b3 = gi13*b_1 + gi23*b_2 + gi33*b_3;

  const Real beta2 =
      b_1*b1 + b_2*b2 + b_3*b3;

  const Real alpha_sq =
      -g(I00,i) + beta2;

  const Real a2i = 1.0/alpha_sq;

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

void Metric::ConstructCellCovariantMetric(
    int k, int j, int i,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const {

  Coordinates *pcoord = pmy_block->pcoord;

  const Real theta = pcoord->x2v(j);
  const Real r = pcoord->x1v(i);

  const Real sintheta = std::sin(theta);
  const Real f = 1.0 - 2.0*bh_mass_/r;

  g00 = -f + 2.0*delta_m_(i)/r + 2.0*f*Psi_(i);
  g01 = 0.0;
  g02 = 0.0;
  g03 = 0.0;
  
  g11 = 1.0/f + 2.0*delta_m_(i)/(r*f*f);
  g12 = 0.0;
  g13 = 0.0;
  
  g22 = r*r;
  g23 = 0.0;
  g33 = r*r*sintheta*sintheta;
}

// C++ headers
#include <cmath>
#include <sstream>
#include <limits>
#include <algorithm>
#include <iomanip>

#include "metric.hpp"

#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#include "../hydro/hydro.hpp"

#include "../gravity/dynamic_metric_driver.hpp"

Metric::Metric(MeshBlock *pmb, ParameterInput *pin)
  : pmy_block(pmb) {
}

Metric::~Metric() {
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
  pmy_block->pmy_mesh->pmetric_driver->CellMetric(
        pmy_block,
	k,j,il,iu,
	g, g_inv);
  for(int i=il; i<=iu; ++i){
    InvertMetric(i, g, g_inv);
  }
  
}

void Metric::Face1Metric(const int k, const int j, const int il, const int iu,
                                AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  pmy_block->pmy_mesh->pmetric_driver->Face1Metric(
        pmy_block,
	k,j,
	il,iu,
	g, g_inv);
}

void Metric::Face2Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  pmy_block->pmy_mesh->pmetric_driver->Face2Metric(
        pmy_block,
	k,j,
	il,iu,
	g, g_inv);

}

void Metric::Face3Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  pmy_block->pmy_mesh->pmetric_driver->Face3Metric(
        pmy_block,
	k,j,
	il,iu,
	g, g_inv);

}

Real Metric::SqrtMinusG(int k, int j, int i) const {

  Real sqrt_minus_g = pmy_block->pmy_mesh->pmetric_driver->SqrtMinusG(
        pmy_block,k,j,i);
  return sqrt_minus_g;
}

Real Metric::CellDensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->CellDensitizationFactor(pmy_block,k,j,i);
}


Real Metric::Face1DensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face1DensitizationFactor(pmy_block,k,j,i);
}


Real Metric::Face2DensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face2DensitizationFactor(pmy_block,k,j,i);
}


Real Metric::Face3DensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face3DensitizationFactor(pmy_block,k,j,i);
}

Real Metric::CellSpatialDensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->CellSpatialDensitizationFactor(pmy_block,k,j,i);
}

Real Metric::Face1SpatialDensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face1SpatialDensitizationFactor(pmy_block,k,j,i);
}


Real Metric::Face2SpatialDensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face2SpatialDensitizationFactor(pmy_block,k,j,i);
}


Real Metric::Face3SpatialDensitizationFactor(int k, int j, int i) const {
  return pmy_block->pmy_mesh->pmetric_driver->Face3SpatialDensitizationFactor(pmy_block,k,j,i);
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


  pmy_block->pmy_mesh->pmetric_driver->ConstructCellCovariantMetric(
        pmy_block,k,j,i,
	g00, g01, g02, g03,
	g11, g12, g13,
	g22, g23, g33);

}

void Metric::CellMetricRadialDerivatives(const int k, const int j,
					 const int il, const int iu,
					 AthenaArray<Real> &d1_g00,
					 AthenaArray<Real> &d1_g11) const {
  pmy_block->pmy_mesh->pmetric_driver->CellMetricRadialDerivatives(
        pmy_block,k,j,il,iu, d1_g00, d1_g11);
}

// C++ headers
#include <cmath>
#include <sstream>
#include <limits>
#include <algorithm>

#include "metric.hpp"

#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#include "../hydro/hydro.hpp"

#include "../gravity/dynamic_metric_driver.hpp"

Metric::Metric(MeshBlock *pmb, ParameterInput *pin)
  : pmy_block(pmb) {

  // const Real bh_mass = GetBlackHoleMass();
  // const Real bh_spin = GetBlackHoleSpin();
  
  // bh_mass_prev_    = bh_mass;
  // bh_mass_pending_ = bh_mass;
  // mdot_bh_         = 0.0;
  
  // bh_spin_prev_    = bh_spin;
  // bh_spin_pending_ = bh_spin;
  // angdot_bh_       = 0.0;
  
}

Metric::~Metric() {
}

// Real Metric::BlackHoleMassAccretionRate(const AthenaArray<Real> &x1flux) const {
  
//   MeshBlock *pmb = pmy_block;
//   Coordinates *pcoord = pmb->pcoord;

//   // This MeshBlock does not touch the physical inner-x1 boundary.
//   if (pmb->pbval->block_bcs[BoundaryFace::inner_x1]
//       == BoundaryFlag::block) {
//     return 0.0;
//   }

//   Real mdot = 0.0;

//   for (int k=pmb->ks; k<=pmb->ke; ++k) {
//     for (int j=pmb->js; j<=pmb->je; ++j) {
//       const Real area = pcoord->GetFace1Area(k, j, pmb->is);
      
//       // inward flux is negative
//       mdot -= area * x1flux(IDN, k, j, pmb->is);
//     }
//   }

//   const Real mass_to_length =
//       pmb->pmy_mesh->punit->grav_const_code
//       / SQR(pmb->pmy_mesh->punit->speed_of_light_code);

//   return mass_to_length * mdot;

// }

// void Metric::CommitBlackHoleMass(){
//   BlackHoleMassStorage() = bh_mass_pending_;
// }


AthenaArray<Real>& Metric::PsiFace1() {
  return pmy_block->ruser_meshblock_data[0];
}

const AthenaArray<Real>& Metric::PsiFace1() const {
  return pmy_block->ruser_meshblock_data[0];
}

AthenaArray<Real>& Metric::DeltaMFace1() {
  return pmy_block->ruser_meshblock_data[1];
}

const AthenaArray<Real>& Metric::DeltaMFace1() const {
  return pmy_block->ruser_meshblock_data[1];
}

Real Metric::CellPsi(int i) const {
  const auto &psi = PsiFace1();
  return 0.5*(psi(i) + psi(i+1));
}

Real Metric::CellDeltaM(int i) const {
  const auto &dm = DeltaMFace1();
  return 0.5*(dm(i) + dm(i+1));
}


// // delta_m_ and Psi_ are derived from fluid distribution.
// void Metric::Update() {

//   auto &Psi_face1 = PsiFace1();
//   auto &delta_m_face1 = DeltaMFace1();
  
//   const Real bh_mass = GetBlackHoleMass();

//   Coordinates *pcoord = pmy_block->pcoord;
//   Hydro *phydro = pmy_block->phydro;

//   const Real mass_to_length =
//     pmy_block->pmy_mesh->punit->grav_const_code
//     / SQR(pmy_block->pmy_mesh->punit->speed_of_light_code);
  
//   const int is = pmy_block->is;
//   const int ie = pmy_block->ie;
//   const int js = pmy_block->js;
//   const int je = pmy_block->je;
//   const int ks = pmy_block->ks;
//   const int ke = pmy_block->ke;
  
//   const int nc1 = pmy_block->ncells1;
//   //const int nc2 = pmy_block->ncells2;
//   //const int nc3 = pmy_block->ncells3;

//   AthenaArray<Real> dm_shell;
//   dm_shell.NewAthenaArray(nc1);
//   dm_shell.ZeroClear();

//   AthenaArray<Real> vol;
//   vol.NewAthenaArray(nc1);

//   for (int k=ks; k<=ke; ++k) {
//     for (int j=js; j<=je; ++j) {
//       pcoord->CellVolume(k, j, is, ie, vol);
// #pragma omp simd
//       for (int i=is; i<=ie; ++i) {
//         dm_shell(i) += phydro->w(IDN, k, j, i) * vol(i) * mass_to_length;
//       }
//     }
//   }

//   std::cout << "pmetric dm_shell:" << std::endl;
//   for(int i=is; i<=ie; ++i){
//     printf("i, dm_shell = %5d %25.16e\n", i, dm_shell(i));
//   }


//   delta_m_face1(is) = 0.0;

//   for (int i=is; i<=ie; ++i) {
//     delta_m_face1(i+1) = delta_m_face1(i) + dm_shell(i);
//   }

//   std::cout << "pmetric delta_m_face1:" << std::endl;
//   for(int i=is; i<=ie; ++i){
//     printf("i, delta_m_face1 = %5d %25.16e\n", i, delta_m_face1(i));
//   }

  
//   dm_shell.ZeroClear();

//   for (int k=ks; k<=ke; ++k) {
//     for (int j=js; j<=je; ++j) {
//       pcoord->CellVolume(k, j, is, ie, vol);
// #pragma omp simd
//       for (int i=is; i<=ie; ++i) {
//         const Real r = pcoord->x1v(i);
//         dm_shell(i) += phydro->w(IDN, k, j, i)/(r-2.0*bh_mass) * vol(i) * mass_to_length;
//       }
//     }
//   }

//   std::cout << "pmetric integrant of Psi:" << std::endl;
//   for(int i=is; i<=ie; ++i){
//     printf("i, integr Psi = %5d %25.16e\n", i, dm_shell(i));
//   }

  
//   Psi_face1(ie+1) = 0.0;

//   for (int i=ie; i>=is; --i) {
//     Psi_face1(i) = Psi_face1(i+1) + dm_shell(i);
//   }

//   std::cout << "pmetric Psi_face1:" << std::endl;
//   for(int i=is; i<=ie; ++i){
//     printf("i, Psi = %5d %25.16e\n", i, Psi_face1(i));
//   }


//   // inner radial ghost faces
//   for (int i=is-1; i>=0; --i) {
//     delta_m_face1(i) = delta_m_face1(is);
//     Psi_face1(i)     = Psi_face1(is);
//   }
  
//   // outer radial ghost faces
//   for (int i=ie+2; i<=nc1; ++i) {
//     delta_m_face1(i) = delta_m_face1(ie+1);
//     Psi_face1(i)     = Psi_face1(ie+1);
//   }
  
// }

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
  
  for(int i=il; i<=iu; ++i){
    const Real r = pcoord->x1v(i);
    
    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;
    ConstructCellCovariantMetric(k,j,i,
      g00, g01, g02, g03,
      g11, g12, g13, g22, g23, g33);
    
    g(I00, i) = g00;
    g(I01, i) = g01;
    g(I02, i) = g02;
    g(I03, i) = g03;
    g(I11, i) = g11;
    g(I12, i) = g12;
    g(I13, i) = g13;
    g(I22, i) = g22;
    g(I23, i) = g23;
    g(I33, i) = g33;
    
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

Real Metric::CellDensitizationFactor(int k, int j, int i) const {
  // return 1.0;
  const Real sqrt_minus_g = SqrtMinusG(k,j,i);
  const Real r = pmy_block->pcoord->x1v(i);
  const Real theta = pmy_block->pcoord->x2v(j);
  return sqrt_minus_g/(r*r*std::sin(theta));
}


Real Metric::Face1DensitizationFactor(int k, int j, int i) const {
  Coordinates *pcoord = pmy_block->pcoord;

  auto &Psi_face1 = PsiFace1();
  auto &delta_m_face1 = DeltaMFace1();
  
  const Real r = pcoord->x1f(i);
  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3v(k);
  const Real Psi = Psi_face1(i);
  const Real dm = delta_m_face1(i);

  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  ConstructCovariantMetric(r, theta, phi, Psi, dm,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);
}


Real Metric::Face2DensitizationFactor(int k, int j, int i) const {
  Coordinates *pcoord = pmy_block->pcoord;

  const Real r = pcoord->x1v(i);
  const Real theta = pcoord->x2f(j);
  const Real phi = pcoord->x3v(k);
  const Real Psi = CellPsi(i);
  const Real dm = CellDeltaM(i);

  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  ConstructCovariantMetric(r, theta, phi, Psi, dm,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);
}


Real Metric::Face3DensitizationFactor(int k, int j, int i) const {
  Coordinates *pcoord = pmy_block->pcoord;

  const Real r = pcoord->x1v(i);
  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3f(k);
  const Real Psi = CellPsi(i);
  const Real dm = CellDeltaM(i);

  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  ConstructCovariantMetric(r, theta, phi, Psi, dm,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);
}

Real Metric::GetBlackHoleMass() const {
  return pmy_block->pmy_mesh->ruser_mesh_data[0](0);
}
Real Metric::GetBlackHoleSpin() const {
  return pmy_block->pmy_mesh->ruser_mesh_data[1](0);
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

  //Coordinates *pcoord = pmy_block->pcoord;

  // const Real r = pcoord->x1v(i);
  // const Real theta = pcoord->x2v(j);
  // const Real phi = pcoord->x3v(k);
  
  // ConstructCovariantMetric(r, theta, phi, CellPsi(i), CellDeltaM(i),
  //     g00, g01, g02, g03,
  //     g11, g12, g13, g22, g23, g33);
  
}

void Metric::ConstructCovariantMetric(
    Real r, Real theta, Real phi, Real Psi, Real dm,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const {

  const Real bh_mass = GetBlackHoleMass();
  gravity_model_.ConstructCovariantMetric(
	r, theta, phi, Psi, dm, bh_mass,
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

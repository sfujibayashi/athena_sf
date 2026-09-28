
// C++ headers
#include <sstream>

// Athena++ headers
#include "monopole_gravity.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../metric/metric.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

namespace{

  int GlobalRadialIndex(const MeshBlock *pmb, int i){
    int ig = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
    return ig;
  }
  
}

MonopoleGravity::MonopoleGravity(Mesh *pm, ParameterInput *pin)
  : pmy_mesh_(pm) {

  if (pm->multilevel) {
    std::stringstream msg;
    msg << "### FATAL ERROR in monopole_gravity.cpp" << std::endl
        << "Multi-level is not supported."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  
  nr_ = pm->mesh_size.nx1;

  dm_shell_global_.NewAthenaArray(nr_);
  delta_m_face_global_.NewAthenaArray(nr_ + 1);
  Psi_face_global_.NewAthenaArray(nr_ + 1);
  
  dm_shell_global_.ZeroClear();
  delta_m_face_global_.ZeroClear();
  Psi_face_global_.ZeroClear();
}

MonopoleGravity::~MonopoleGravity() {
  dm_shell_global_.DeleteAthenaArray();
  delta_m_face_global_.DeleteAthenaArray();
  Psi_face_global_.DeleteAthenaArray();
}


void MonopoleGravity::Update(){

  const Real bh_mass = GetBlackHoleMass();
  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmy_mesh_->my_blocks(0)->ncells1);
  
  dm_shell_global_.ZeroClear();
  
  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmy_mesh_->my_blocks(b)->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          
          int ig = GlobalRadialIndex(pmb, i);

          dm_shell_global_(ig) +=
            pmb->phydro->w(IDN,k,j,i)
            * vol(i);
        }
      }
    }
  }
  
  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  for(int i=0; i<nr_; ++i){
    dm_shell_global_(i) *= mass_to_length;
  }
  
  std::cout << "pmonotgrav dm_shell:" << std::endl;
  for(int i=0; i<nr_; ++i){
    printf("ig, dm_shell = %5d %25.16e\n", i, dm_shell_global_(i));
  }
  
  delta_m_face_global_(0) = 0.0;
  for (int i=0; i<nr_; ++i) {
    delta_m_face_global_(i+1) = delta_m_face_global_(i) + dm_shell_global_(i);
  }

  std::cout << "pmonograv delta_m_face1:" << std::endl;
  for(int i=0; i<=nr_; ++i){
    printf("i, delta_m_face1 = %5d %25.16e\n", i, delta_m_face_global_(i));
  }
  
  dm_shell_global_.ZeroClear();
  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmy_mesh_->my_blocks(b)->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          
          int ig = GlobalRadialIndex(pmb, i);

          const Real r = pmb->pcoord->x1v(i);
          dm_shell_global_(i) += 
            pmb->phydro->w(IDN,k,j,i)/(r-2.0*bh_mass) * vol(i) * mass_to_length;
        }
      }
    }
  }

  std::cout << "pmonograv integrant of Psi:" << std::endl;
  for(int i=0; i<nr_; ++i){
    printf("i, integr Psi = %5d %25.16e\n", i, dm_shell_global_(i));
  }
  
  Psi_face_global_(nr_) = 0.0;

  for (int i=nr_; i>=0; --i) {
    Psi_face_global_(i) = Psi_face_global_(i+1) + dm_shell_global_(i);
  }

  std::cout << "pmonotgrav Psi_face1:" << std::endl;
  for(int i=0; i<=nr_; ++i){
    printf("i, Psi = %5d %25.16e\n", i, Psi_face_global_(i));
  }
}


Real& MonopoleGravity::BlackHoleMassStorage() {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

const Real& MonopoleGravity::BlackHoleMassStorage() const {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

Real& MonopoleGravity::BlackHoleSpinStorage() {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

const Real& MonopoleGravity::BlackHoleSpinStorage() const {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

Real MonopoleGravity::GetBlackHoleMass() const {
  return BlackHoleMassStorage();
}
Real MonopoleGravity::GetBlackHoleSpin() const {
  return BlackHoleSpinStorage();
}

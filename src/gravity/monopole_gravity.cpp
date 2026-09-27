#include "monopole_gravity.hpp"

#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../metric/metric.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

MonopoleGravity::MonopoleGravity(Mesh *pm, ParameterInput *pin)
  : pmy_mesh_(pm) {

  nr_ = pm->mesh_size.nx1;

  dm_shell_global_.NewAthenaArray(nr_);
  delta_m_face_global_.NewAthenaArray(nr_ + 1);
  psi_face_global_.NewAthenaArray(nr_ + 1);
  
  dm_shell_global_.ZeroClear();
  delta_m_face_global_.ZeroClear();
  psi_face_global_.ZeroClear();
}

MonopoleGravity::~MonopoleGravity() {
  dm_shell_global_.DeleteAthenaArray();
  delta_m_face_global_.DeleteAthenaArray();
  psi_face_global_.DeleteAthenaArray();
}

int GlobalRadialIndex(MeshBlock *pmb, int i){
  int ig = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
  return ig;
}

void MonopoleGravity::Update(){

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
}


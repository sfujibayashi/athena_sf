
// C++ headers
#include <sstream>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../metric/metric.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

#include "monopole_gravity_driver.hpp"
#include "dynamic_metric_driver.hpp"

namespace{

  int GlobalRadialIndex(const MeshBlock *pmb, int i){
    int ig = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
    return ig;
  }
  
}

MonopoleGravityDriver::MonopoleGravityDriver(Mesh *pm, ParameterInput *pin)
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

  
  const Real bh_mass = GetBlackHoleMass();
  const Real bh_spin = GetBlackHoleSpin();
  
  bh_mass_prev_    = bh_mass;
  bh_mass_pending_ = bh_mass;
  mdot_bh_         = 0.0;
  
  bh_spin_prev_    = bh_spin;
  bh_spin_pending_ = bh_spin;
  angdot_bh_       = 0.0;

}

MonopoleGravityDriver::~MonopoleGravityDriver() {
  dm_shell_global_.DeleteAthenaArray();
  delta_m_face_global_.DeleteAthenaArray();
  Psi_face_global_.DeleteAthenaArray();
}


void MonopoleGravityDriver::Update(){

  const Real bh_mass = GetBlackHoleMass();

  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmy_mesh_->block_size.nx1+2*NGHOST);
  
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

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
                dm_shell_global_.data(),
                nr_,
                MPI_ATHENA_REAL,
                MPI_SUM,
                MPI_COMM_WORLD);
#endif
  
  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  for(int i=0; i<nr_; ++i){
    dm_shell_global_(i) *= mass_to_length;
  }
  
  delta_m_face_global_(0) = 0.0;
  for (int i=0; i<nr_; ++i) {
    delta_m_face_global_(i+1) = delta_m_face_global_(i) + dm_shell_global_(i);
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
          dm_shell_global_(ig) += 
            pmb->phydro->w(IDN,k,j,i)/(r-2.0*bh_mass) * vol(i) * mass_to_length;
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
                dm_shell_global_.data(),
                nr_,
                MPI_ATHENA_REAL,
                MPI_SUM,
                MPI_COMM_WORLD);
#endif
  
  Psi_face_global_(nr_) = 0.0;
  for (int i=nr_-1; i>=0; --i) {
    Psi_face_global_(i) = Psi_face_global_(i+1) + dm_shell_global_(i);
  }

  // Provide global delta_m and Psi to MeshBlock-local storage
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    auto &dm  = pmb->pmetric->DeltaMFace1();
    auto &Psi = pmb->pmetric->PsiFace1();
    
    for (int i=0; i<=pmb->ncells1; ++i) {
      int igf = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
      
      if (igf < 0) {
        dm(i)  = delta_m_face_global_(0);
        Psi(i) = Psi_face_global_(0);
      } else if (igf > nr_) {
        dm(i)  = delta_m_face_global_(nr_);
        Psi(i) = Psi_face_global_(nr_);
      } else {
        dm(i)  = delta_m_face_global_(igf);
        Psi(i) = Psi_face_global_(igf);
      }
    }
  }
  
}


Real& MonopoleGravityDriver::BlackHoleMassStorage() {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

const Real& MonopoleGravityDriver::BlackHoleMassStorage() const {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

Real& MonopoleGravityDriver::BlackHoleSpinStorage() {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

const Real& MonopoleGravityDriver::BlackHoleSpinStorage() const {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

Real MonopoleGravityDriver::GetBlackHoleMass() const {
  return BlackHoleMassStorage();
}
Real MonopoleGravityDriver::GetBlackHoleSpin() const {
  return BlackHoleSpinStorage();
}

Real MonopoleGravityDriver::BlackHoleMassAccretionRate() const {
  
  Real mdot = 0.0;
  
  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    Coordinates *pcoord = pmb->pcoord;

    // This MeshBlock does not touch the physical inner-x1 boundary.
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::block) {
      continue;
    }

    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        const Real area = pcoord->GetFace1Area(k, j, pmb->is);
        
        // inward flux is negative
        mdot -= area * pmb->phydro->flux[X1DIR](IDN, k, j, pmb->is);
      }
    }
  }
  
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
                &mdot,
                1,
                MPI_ATHENA_REAL,
                MPI_SUM,
                MPI_COMM_WORLD);
#endif
  
  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  return mass_to_length * mdot;
  
}


void MonopoleGravityDriver::UpdateBlackHoleMass(int stage){
   
  mdot_bh_ = BlackHoleMassAccretionRate();

  if (stage == 1) {
    bh_mass_prev_ = GetBlackHoleMass();
    bh_mass_pending_ = bh_mass_prev_ + 0.5*pmy_mesh_->dt*mdot_bh_;
  }

  if (stage == 2) {
    bh_mass_pending_ = bh_mass_prev_ + pmy_mesh_->dt*mdot_bh_;
  }

  BlackHoleMassStorage() = bh_mass_pending_;
}

Real MonopoleGravityDriver::GetBlackHoleMassAccretionRate() const {
  return mdot_bh_;
}

void MonopoleGravityDriver::CellMetricRadialDerivatives(
         MeshBlock *pmb,
         const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &d1_g00,
	 AthenaArray<Real> &d1_g11) const{
  
  const auto &Psi_face1 = PsiFace1(pmb);
  const auto &delta_m_face1 = DeltaMFace1(pmb);
  const Real bh_mass = GetBlackHoleMass();
  const Real theta = pmb->pcoord->x2v(j);
  const Real phi = pmb->pcoord->x3v(k);

  for (int i = il; i <= iu; ++i) {

    const Real r = pmb->pcoord->x1v(i);
    const Real Psi = CellPsi(pmb,i);
    const Real delta_m = CellDeltaM(pmb,i);
    const Real dxf = pmb->pcoord->x1f(i+1)-pmb->pcoord->x1f(i);
    const Real d1_Psi = (Psi_face1(i+1) - Psi_face1(i))/dxf;
    const Real d1_delta_m = (delta_m_face1(i+1) - delta_m_face1(i))/dxf;
    
    gravity_model_.MetricRadialDerivatives(
       r, theta, phi, Psi, delta_m, bh_mass,
       d1_Psi, d1_delta_m,
       d1_g00(i), d1_g11(i));

  }

}  

AthenaArray<Real>& MonopoleGravityDriver::PsiFace1(MeshBlock *pmb) {
  return pmb->ruser_meshblock_data[0];
}

const AthenaArray<Real>& MonopoleGravityDriver::PsiFace1(MeshBlock *pmb) const {
  return pmb->ruser_meshblock_data[0];
}

AthenaArray<Real>& MonopoleGravityDriver::DeltaMFace1(MeshBlock *pmb) {
  return pmb->ruser_meshblock_data[1];
}

const AthenaArray<Real>& MonopoleGravityDriver::DeltaMFace1(MeshBlock *pmb) const {
  return pmb->ruser_meshblock_data[1];
}

Real MonopoleGravityDriver::CellPsi(MeshBlock *pmb, int i) const {
  const auto &psi = PsiFace1(pmb);
  return 0.5*(psi(i) + psi(i+1));
}

Real MonopoleGravityDriver::CellDeltaM(MeshBlock *pmb, int i) const {
  const auto &dm = DeltaMFace1(pmb);
  return 0.5*(dm(i) + dm(i+1));
}


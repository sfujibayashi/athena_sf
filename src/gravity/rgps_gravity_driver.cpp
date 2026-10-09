
// C++ headers
#include <sstream>
#include <limits>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../metric/metric.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../eos/eos.hpp"

#include "rgps_gravity.hpp"
#include "rgps_gravity_driver.hpp"
#include "dynamic_metric_driver.hpp"

namespace{

  int GlobalRadialIndex(const MeshBlock *pmb, int i){
    int ig = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
    return ig;
  }
  
}

RGPSGravityDriver::RGPSGravityDriver(Mesh *pm, ParameterInput *pin)
  : pmy_mesh_(pm) {

  if (pm->multilevel) {
    std::stringstream msg;
    msg << "### FATAL ERROR in RGPSGravityDriver::RGPSGravityDriver" << std::endl
        << "Multi-level is not supported."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  if (pmy_mesh_->nblocal != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in RGPSGravityDriver::RGPSGravityDriver"
        << std::endl
        << "RGPS dynamic metric currently requires exactly one "
        << "MeshBlock per MPI rank." << std::endl
        << "nblocal = " << pmy_mesh_->nblocal << std::endl;
    ATHENA_ERROR(msg);
  }
  
  nr_ = pm->mesh_size.nx1;

  r_cell_global_.NewAthenaArray(nr_);
  r_face_global_.NewAthenaArray(nr_+1);
  vol_.NewAthenaArray(pmy_mesh_->block_size.nx1+2*NGHOST);

  InitializeRadialGrid();

  calE_shell_global_.NewAthenaArray(nr_);
  Srr_shell_global_.NewAthenaArray(nr_);
  // Sr_shell_global_.NewAthenaArray(nr_);
  dmgrav_dr_cell_global_.NewAthenaArray(nr_);
  // dmgrav_dt_global_.NewAthenaArray(nr_);
  dPhi_dr_cell_global_.NewAthenaArray(nr_);
  X_sq_cell_global_.NewAthenaArray(nr_);
  mgrav_face_global_.NewAthenaArray(nr_ + 1);
  Phi_face_global_.NewAthenaArray(nr_ + 1);

  calE_shell_global_.ZeroClear();
  Srr_shell_global_.ZeroClear();
  // Sr_shell_global_.ZeroClear();
  dmgrav_dr_cell_global_.ZeroClear();
  // dmgrav_dt_global_.ZeroClear();
  dPhi_dr_cell_global_.ZeroClear();
  X_sq_cell_global_.ZeroClear();
  mgrav_face_global_.ZeroClear();
  Phi_face_global_.ZeroClear();
  
  const Real bh_mass = GetBlackHoleMass();
  const Real bh_spin = GetBlackHoleSpin();
  
  bh_mass_prev_    = bh_mass;
  bh_mass_pending_ = bh_mass;
  mdot_bh_         = 0.0;
  
  bh_spin_prev_    = bh_spin;
  bh_spin_pending_ = bh_spin;
  angdot_bh_       = 0.0;

}

RGPSGravityDriver::~RGPSGravityDriver() {
  calE_shell_global_.DeleteAthenaArray();
  Srr_shell_global_.DeleteAthenaArray();
  // Sr_shell_global_.DeleteAthenaArray();
  dmgrav_dr_cell_global_.DeleteAthenaArray();
  // dmgrav_dt_global_.DeleteAthenaArray();
  dPhi_dr_cell_global_.DeleteAthenaArray();
  mgrav_face_global_.DeleteAthenaArray();
  Phi_face_global_.DeleteAthenaArray();
  X_sq_cell_global_.DeleteAthenaArray();
}

void RGPSGravityDriver::InitializeRadialGrid(){

  AthenaArray<Real> r_sum, rf_sum;
  AthenaArray<int> count, countf;

  r_sum.NewAthenaArray(nr_);
  rf_sum.NewAthenaArray(nr_+1);
  count.NewAthenaArray(nr_);
  countf.NewAthenaArray(nr_+1);

  r_sum.ZeroClear();
  rf_sum.ZeroClear();
  count.ZeroClear();
  countf.ZeroClear();
  
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      const int ig = GlobalRadialIndex(pmb, i);
      r_sum(ig) += pmb->pcoord->x1v(i);
      count(ig)++;
    }
    
    for (int i=pmb->is; i<=pmb->ie+1; ++i) {
      const int igf = pmb->loc.lx1*pmb->block_size.nx1 + (i-pmb->is);
      rf_sum(igf) += pmb->pcoord->x1f(i);
      countf(igf)++;
    }
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, r_sum.data(), nr_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, rf_sum.data(), nr_+1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count.data(), nr_, MPI_ATHENA_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, countf.data(), nr_+1, MPI_ATHENA_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  for(int i=0; i<nr_;++i){
    r_cell_global_(i) = r_sum(i)/count(i);
  }

  for(int i=0; i<nr_+1;++i){
    r_face_global_(i) = rf_sum(i)/countf(i);
  }
}

void RGPSGravityDriver::ConstructMgravFromPrimitive(){

  const Real bh_mass = GetBlackHoleMass();
  
  AthenaArray<Real> Em_shell_global;
  AthenaArray<Real> Er_shell_global;
  AthenaArray<Real> Ea_shell_global;

  Em_shell_global.NewAthenaArray(nr_);
  Er_shell_global.NewAthenaArray(nr_);
  Ea_shell_global.NewAthenaArray(nr_);

  Em_shell_global.ZeroClear();
  Er_shell_global.ZeroClear();
  Ea_shell_global.ZeroClear();

  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmy_mesh_->my_blocks(b)->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
	const Real theta = pmy_mesh_->my_blocks(b)->pcoord->x2v(j);
	const Real sintheta = std::sin(theta);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
	  const Real r = pmy_mesh_->my_blocks(b)->pcoord->x1v(i);

	  const Real rho  = pmb->phydro->w(IDN,k,j,i);
	  const Real pgas = pmb->phydro->w(IPR,k,j,i);
	  const Real uu1  = pmb->phydro->w(IVX,k,j,i);
	  const Real uu2  = pmb->phydro->w(IVY,k,j,i);
	  const Real uu3  = pmb->phydro->w(IVZ,k,j,i);

	  const Real egas = pmb->peos->EgasFromRhoP(rho, pgas);
	  const Real rhoh = rho + egas + pgas;
	  
	  const Real kr = rhoh*uu1*uu1;
	  const Real ka = rhoh*r*r*(uu2*uu2 + uu3*uu3*sintheta*sintheta);
	  
          int ig = GlobalRadialIndex(pmb, i);
          Em_shell_global(ig) += (rhoh-pgas) * vol_(i);
          Er_shell_global(ig) += kr * vol_(i);
          Ea_shell_global(ig) += ka * vol_(i);
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, Em_shell_global.data(), nr_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, Er_shell_global.data(), nr_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, Ea_shell_global.data(), nr_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  // Get angular average value by dividing with the shell volume.
  for(int i=0; i<nr_; ++i){
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);

    const Real dV = 4.0*M_PI/3.0 * (rp*rp*rp - rm*rm*rm);
    Em_shell_global(i) /= dV;
    Er_shell_global(i) /= dV;
    Ea_shell_global(i) /= dV;
  }
  
  // E has code units of mass density
  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  mgrav_face_global_(0) = bh_mass;
  for (int i=0; i<nr_; ++i) {
    const Real mL = mgrav_face_global_(i);
    const Real rc = r_cell_global_(i);
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);

    const Real dVm= 4.0*M_PI/3.0 * (rc*rc*rc - rm*rm*rm);
    
    const Real Eatmos = pmy_mesh_->my_blocks(0)->peos->GetEnergyFloor(r_cell_global_(i));

    // Energy in units of code length.
    const Real dEr_len    = Er_shell_global(i)*dVm * mass_to_length;
    const Real dErest_len = (Em_shell_global(i) + Ea_shell_global(i) - Eatmos)*dVm * mass_to_length;
    
    const Real b = 1.0 - 2.0/rc*(mL + dErest_len);
    const Real c = 2.0*dEr_len/rc;

    const Real disc = b*b - 4.0*c;
    if (disc < 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in RGPSGravityDriver::ConstructMgravFromPrimitive" << std::endl
	  << "b^2-4c < 0."
	  << std::endl;
      ATHENA_ERROR(msg);
    }
    
    const Real Xinv_sq = 0.5*(b + std::sqrt(disc));
    
    if (Xinv_sq <= 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in RGPSGravityDriver::ConstructMgravFromPrimitive" << std::endl
	  << "Xinv_sq <= 0."
	  << std::endl;
      ATHENA_ERROR(msg);
    }

    const Real X_sq = 1.0/Xinv_sq;
    X_sq_cell_global_(i) = X_sq;

    const Real dV = 4.0*M_PI/3.0 * (rp*rp*rp - rm*rm*rm);
    const Real Egrav = Em_shell_global(i) + X_sq*Er_shell_global(i) + Ea_shell_global(i) - Eatmos;
    const Real dm_full = dV*Egrav*mass_to_length;
    mgrav_face_global_(i+1) = mL + dm_full;
    // mgrav_face_global_(i+1) = mL + std::max(dcalE_len*Xinv - dEat_len, 0.0);

    dmgrav_dr_cell_global_(i) = 4.0*M_PI*rc*rc*Egrav*mass_to_length;
  }
  

  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    auto &mgrav_face1=MgravFace1(pmb);

    for (int i=0; i<=pmb->ncells1; ++i) {
      int igf = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
      
      if (igf < 0) {
        mgrav_face1(i)  = mgrav_face_global_(0);
      } else if (igf > nr_) {
        mgrav_face1(i)  = mgrav_face_global_(nr_);
      } else {
        mgrav_face1(i)  = mgrav_face_global_(igf);
      }
    }
  } 
}

void RGPSGravityDriver::ConstructPhiFromPrimitive(){

  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  AthenaArray<Real> Sr_r_shell_global;
  Sr_r_shell_global.NewAthenaArray(nr_);
  Sr_r_shell_global.ZeroClear();
  
  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmy_mesh_->my_blocks(b)->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
        for (int i=pmb->is; i<=pmb->ie; ++i) {

	  const Real r = pmb->pcoord->x1v(i);

	  const Real rho  = pmb->phydro->w(IDN,k,j,i);
	  const Real pgas = pmb->phydro->w(IPR,k,j,i);
	  const Real uu1  = pmb->phydro->w(IVX,k,j,i);
	  const Real uu2  = pmb->phydro->w(IVY,k,j,i);
	  const Real uu3  = pmb->phydro->w(IVZ,k,j,i);

	  const int ig = GlobalRadialIndex(pmb, i);
	  const Real X_sq  = X_sq_cell_global_(ig);

	  const Real egas = pmb->peos->EgasFromRhoP(rho, pgas);
	  const Real wtot = rho + egas + pgas; // rho*h
	  const Real Sr_r = wtot*X_sq*uu1*uu1 + pgas;

          Sr_r_shell_global(ig) += Sr_r*vol_(i);
        }
      }
    }
  }
  
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, Sr_r_shell_global.data(), nr_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  // Get angular average value by dividing with the shell volume.
  for(int i=0; i<nr_; ++i){
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);

    const Real dV = 4.0*M_PI/3.0 * (rp*rp*rp - rm*rm*rm);
    Sr_r_shell_global(i) /= dV;
  }
  

  const Real mout = mgrav_face_global_(nr_);
  const Real rout = r_face_global_(nr_);
  Phi_face_global_(nr_) = 0.5*std::log(1.0 - 2.0*mout/rout);
  for (int i=nr_-1; i>=0; --i) {
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);
    const Real dr = rp-rm;
    
    // cell-centered values
    const Real r = r_cell_global_(i);
    const Real X_sq = X_sq_cell_global_(i);
    const Real mgrav = 0.5*r*(1.0-1.0/X_sq);
    const Real Sr_r = Sr_r_shell_global(i);

    const Real Patm = pmy_mesh_->my_blocks(0)->peos->GetPressureFloor(r);
    
    const Real dPhi_dr = X_sq*(mgrav/(r*r) + 4.0*M_PI*r*(Sr_r-Patm)*mass_to_length);
    Phi_face_global_(i) = Phi_face_global_(i+1) - dr * dPhi_dr;
    
    dPhi_dr_cell_global_(i) = dPhi_dr;
  }

  
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    auto &Phi_face1=PhiFace1(pmb);

    for (int i=0; i<=pmb->ncells1; ++i) {
      int igf = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
      
      if (igf < 0) {
        Phi_face1(i)  = Phi_face_global_(0);
      } else if (igf > nr_) {
        Phi_face1(i)  = Phi_face_global_(nr_);
      } else {
        Phi_face1(i)  = Phi_face_global_(igf);
      }
    }
  } 

}

void RGPSGravityDriver::ConstructMgravFromConserved(){
  
  const Real bh_mass = GetBlackHoleMass();

  // AthenaArray<Real> vol;
  // vol.NewAthenaArray(pmy_mesh_->block_size.nx1+2*NGHOST);
  
  calE_shell_global_.ZeroClear();
  
  // run over MeshBlocks
  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmy_mesh_->my_blocks(b)->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          
          int ig = GlobalRadialIndex(pmb, i);

          calE_shell_global_(ig) += (pmb->phydro->u(IEN,k,j,i)+pmb->phydro->u(IDN,k,j,i)) * vol_(i);
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
                calE_shell_global_.data(),
                nr_,
                MPI_ATHENA_REAL,
                MPI_SUM,
                MPI_COMM_WORLD);
#endif

  // Get angular average value by dividing with the shell volume.
  for(int i=0; i<nr_; ++i){
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);

    const Real dV = 4.0*M_PI/3.0 * (rp*rp*rp - rm*rm*rm);
    calE_shell_global_(i) /= dV;
  }
  
  // calE has code units of mass density
  const Real mass_to_length =
    pmy_mesh_->punit->grav_const_code
    / SQR(pmy_mesh_->punit->speed_of_light_code);
  
  mgrav_face_global_(0) = bh_mass;
  for (int i=0; i<nr_; ++i) {
    const Real mL = mgrav_face_global_(i);
    const Real rc = r_cell_global_(i);
    const Real rm = r_face_global_(i);
    const Real rp = r_face_global_(i+1);

    const Real dVm= 4.0*M_PI/3.0 * (rc*rc*rc - rm*rm*rm);

    const Real Eatmos = pmy_mesh_->my_blocks(0)->peos->GetEnergyFloor(r_cell_global_(i));

    // Energy in units of code length.
    const Real dcalE_len = calE_shell_global_(i)*dVm * mass_to_length;
    const Real dEat_len = Eatmos*dVm * mass_to_length;
    
    const Real b  = dcalE_len/rc;
    const Real cm = 1.0 - 2.0*(mL-dEat_len)/rc;
    
    const Real Xinv = cm/(b + std::sqrt(b*b + cm));
    const Real X_sq = 1.0/(Xinv*Xinv);
    X_sq_cell_global_(i) = X_sq;
    
    const Real dV = 4.0*M_PI/3.0 * (rp*rp*rp - rm*rm*rm);
    const Real Egrav = calE_shell_global_(i)*Xinv - Eatmos;
    const Real dm_full = dV*Egrav*mass_to_length;

    mgrav_face_global_(i+1) = mL + dm_full;
    // mgrav_face_global_(i+1) = mL + std::max(dcalE_len*Xinv - dEat_len, 0.0);

    dmgrav_dr_cell_global_(i) = 4.0*M_PI*rc*rc*Egrav*mass_to_length;
  }

  for (int b=0; b<pmy_mesh_->nblocal; ++b) {
    MeshBlock *pmb = pmy_mesh_->my_blocks(b);
    
    auto &mgrav_face1=MgravFace1(pmb);

    for (int i=0; i<=pmb->ncells1; ++i) {
      int igf = pmb->loc.lx1 * pmb->block_size.nx1 + (i - pmb->is);
      
      if (igf < 0) {
        mgrav_face1(i)  = mgrav_face_global_(0);
      } else if (igf > nr_) {
        mgrav_face1(i)  = mgrav_face_global_(nr_);
      } else {
        mgrav_face1(i)  = mgrav_face_global_(igf);
      }
    }
  } 
}


void RGPSGravityDriver::UpdateBeforeCons2Prim(){
}

void RGPSGravityDriver::UpdateAfterCons2Prim(){
  
}


Real& RGPSGravityDriver::BlackHoleMassStorage() {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

const Real& RGPSGravityDriver::BlackHoleMassStorage() const {
  return pmy_mesh_->ruser_mesh_data[0](0);
}

Real& RGPSGravityDriver::BlackHoleSpinStorage() {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

const Real& RGPSGravityDriver::BlackHoleSpinStorage() const {
  return pmy_mesh_->ruser_mesh_data[1](0);
}

Real RGPSGravityDriver::GetBlackHoleMass() const {
  return BlackHoleMassStorage();
}
Real RGPSGravityDriver::GetBlackHoleSpin() const {
  return BlackHoleSpinStorage();
}

Real RGPSGravityDriver::BlackHoleMassAccretionRate() const {
  
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


void RGPSGravityDriver::UpdateBlackHoleMass(int stage){
   
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

Real RGPSGravityDriver::GetBlackHoleMassAccretionRate() const {
  return mdot_bh_;
}

void RGPSGravityDriver::CellMetricRadialDerivatives(
         MeshBlock *pmb,
         const int k, const int j,
	 const int il, const int iu,
	 AthenaArray<Real> &d1_g00,
	 AthenaArray<Real> &d1_g11) const{
  
  const auto &Phi_face1 = PhiFace1(pmb);
  const auto &mgrav_face1 = MgravFace1(pmb);
  const Real theta = pmb->pcoord->x2v(j);
  const Real phi = pmb->pcoord->x3v(k);

  for (int i = il; i <= iu; ++i) {

    const Real r = pmb->pcoord->x1v(i);
    const Real Phi = CellPhi(pmb,i);

    const int ig = GlobalRadialIndex(pmb, i);
    const Real X_sq = X_sq_cell_global_(ig);
    const Real mgrav = 0.5*r*(1.0-1.0/X_sq);
    
    const Real d1_Phi = dPhi_dr_cell_global_(ig);
    const Real d1_mgrav = dmgrav_dr_cell_global_(ig);
    
    gravity_model_.MetricRadialDerivatives(
       r, theta, phi, Phi, mgrav,
       d1_Phi, d1_mgrav,
       d1_g00(i), d1_g11(i));

  }
}  

void RGPSGravityDriver::Face1Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const {
    // Extract geometric quantities that do not depend on r
  Coordinates *pcoord = pmb->pcoord;

  const auto &Phi_face1 = PhiFace1(pmb);
  const auto &mgrav_face1 = MgravFace1(pmb);

  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3v(k);
  
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

    const Real r = pcoord->x1f(i);
    
    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;

    Real Phi = Phi_face1(i);
    Real mgrav = mgrav_face1(i);

    gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
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
    
    g_inv(I00,i) = 1.0/g00;
    g_inv(I01,i) = 0.0;
    g_inv(I02,i) = 0.0;
    g_inv(I03,i) = 0.0;
    g_inv(I11,i) = 1.0/g11;
    g_inv(I12,i) = 0.0;
    g_inv(I13,i) = 0.0;
    g_inv(I22,i) = 1.0/g22;
    g_inv(I23,i) = 0.0;
    g_inv(I33,i) = 1.0/g33;
  }
  return;

}

void RGPSGravityDriver::Face2Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const {

  Coordinates *pcoord = pmb->pcoord;
  
  const Real theta = pcoord->x2f(j);
  const Real phi = pcoord->x3v(k);
  const bool pole = pcoord->IsPole(j);
  
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

    const Real r = pcoord->x1v(i);

    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;
    
    Real Phi = CellPhi(pmb, i);
    Real mgrav = CellMgrav(pmb, i);
    
    gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
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

    g_inv(I00,i) = 1.0/g00;
    g_inv(I01,i) = 0.0;
    g_inv(I02,i) = 0.0;
    g_inv(I03,i) = 0.0;
    
    g_inv(I11,i) = 1.0/g11;
    g_inv(I12,i) = 0.0;
    g_inv(I13,i) = 0.0;
    
    g_inv(I22,i) = 1.0/g22;
    g_inv(I23,i) = 0.0;

    // Coordinate singularity at theta = 0 or pi.
    // Current Schwarzschild + l=0 self-gravity metric is diagonal.
    if (!pole) {
      g_inv(I33,i) = 1.0/g33;
    }else{
      // g^{phi phi} is singular in spherical coordinates at the pole.
      g_inv(I33,i) = std::numeric_limits<Real>::infinity();
    }
    
  }
}


void RGPSGravityDriver::Face3Metric(
         MeshBlock *pmb,
	 const int k, const int j, 
	 const int il, const int iu,
	 AthenaArray<Real> &g, 
	 AthenaArray<Real> &g_inv) const {

  // Extract geometric quantities that do not depend on r
  Coordinates *pcoord = pmb->pcoord;
  
  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3f(k);
  
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {

    const Real r = pcoord->x1v(i);

    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;
    
    Real Phi = CellPhi(pmb, i);
    Real mgrav = CellMgrav(pmb, i);
    
    gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
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

    g_inv(I00,i) = 1.0/g00;
    g_inv(I01,i) = 0.0;
    g_inv(I02,i) = 0.0;
    g_inv(I03,i) = 0.0;
    g_inv(I11,i) = 1.0/g11;
    g_inv(I12,i) = 0.0;
    g_inv(I13,i) = 0.0;
    g_inv(I22,i) = 1.0/g22;
    g_inv(I23,i) = 0.0;
    g_inv(I33,i) = 1.0/g33;

  }
}


void RGPSGravityDriver::CellMetric(
    MeshBlock *pmb,
    const int k, const int j,
    const int il, const int iu,
    AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv) const {
  
  Coordinates *pcoord = pmb->pcoord;
  
  for(int i=il; i<=iu; ++i){
    const Real r = pcoord->x1v(i);
    
    Real g00, g01, g02, g03;
    Real g11, g12, g13, g22, g23, g33;
    ConstructCellCovariantMetric(
      pmb,
      k,j,i,
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

    g_inv(I00,i) = 1.0/g00;
    g_inv(I01,i) = 0.0;
    g_inv(I02,i) = 0.0;
    g_inv(I03,i) = 0.0;
    g_inv(I11,i) = 1.0/g11;
    g_inv(I12,i) = 0.0;
    g_inv(I13,i) = 0.0;
    g_inv(I22,i) = 1.0/g22;
    g_inv(I23,i) = 0.0;
    g_inv(I33,i) = 1.0/g33;
  }

}

void RGPSGravityDriver::ConstructCellCovariantMetric(
         MeshBlock *pmb,
	 int k, int j, int i,
	 Real &g00, Real &g01, Real &g02, Real &g03,
	 Real &g11, Real &g12, Real &g13,
	 Real &g22, Real &g23, Real &g33) const {

  const Real Phi = CellPhi(pmb, i);
  const Real mgrav  = CellMgrav(pmb, i);

  const Real r = pmb->pcoord->x1v(i);
  const Real theta = pmb->pcoord->x2v(j);
  const Real phi = pmb->pcoord->x3v(k);
  
  gravity_model_.ConstructCovariantMetric(
      r, theta, phi, Phi, mgrav,
      g00, g01, g02, g03,
      g11, g12, g13,
      g22, g23, g33);
    
}

Real RGPSGravityDriver::SqrtMinusG(
         MeshBlock *pmb,
	 const int k, const int j, const int i) const {
  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;

  ConstructCellCovariantMetric(
       pmb, k,j,i,
       g00, g01, g02, g03,
       g11, g12, g13,
       g22, g23, g33);
  
  // Real detgamma = g11*(g22*g33);
  // Real alpha_sq = -g00;
  
  
  return std::sqrt(-g00*(g11*(g22*g33)));
  //return std::sqrt(alpha_sq*detgamma);
}



Real RGPSGravityDriver::CellDensitizationFactor(
        MeshBlock *pmb, int k, int j, int i) const {
  const Real sqrt_minus_g = SqrtMinusG(pmb,k,j,i);
  const Real r = pmb->pcoord->x1v(i);
  const Real theta = pmb->pcoord->x2v(j);
  return sqrt_minus_g/(r*r*std::sin(theta));

  // Real g00, g01, g02, g03;
  // Real g11, g12, g13, g22, g23, g33;

  // ConstructCellCovariantMetric(
  //      pmb, k,j,i,
  //      g00, g01, g02, g03,
  //      g11, g12, g13,
  //      g22, g23, g33);

  //return g11;
}

Real RGPSGravityDriver::Face1DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const {
  Coordinates *pcoord = pmb->pcoord;

  const auto &Phi_face1 = PhiFace1(pmb);
  const auto &mgrav_face1 = MgravFace1(pmb);
  
  const Real r = pcoord->x1f(i);
  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3v(k);
  const Real Phi = Phi_face1(i);
  const Real mgrav = mgrav_face1(i);
  
  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);

}

Real RGPSGravityDriver::Face2DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const {
  Coordinates *pcoord = pmb->pcoord;

  const Real r = pcoord->x1v(i);
  const Real theta = pcoord->x2f(j);
  const Real phi = pcoord->x3v(k);
  const Real Phi = CellPhi(pmb,i);
  const Real mgrav = CellMgrav(pmb,i);
  
  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);

}

Real RGPSGravityDriver::Face3DensitizationFactor(
	MeshBlock *pmb, int k, int j, int i) const {
  Coordinates *pcoord = pmb->pcoord;

  const Real r = pcoord->x1v(i);
  const Real theta = pcoord->x2v(j);
  const Real phi = pcoord->x3f(k);
  const Real Phi = CellPhi(pmb,i);
  const Real mgrav = CellMgrav(pmb,i);
  
  Real g00, g01, g02, g03;
  Real g11, g12, g13, g22, g23, g33;
  
  gravity_model_.ConstructCovariantMetric(r, theta, phi, Phi, mgrav,
        g00, g01, g02, g03,
        g11, g12, g13, g22, g23, g33);

  return std::sqrt((-g00)*g11);

}

Real RGPSGravityDriver::GetEnclosedMassAtInnerBoundary(
    MeshBlock *pmb) const {
  return MgravFace1(pmb)(pmb->is);
}



AthenaArray<Real>& RGPSGravityDriver::PhiFace1(MeshBlock *pmb) {
  return pmb->ruser_meshblock_data[0];
}

const AthenaArray<Real>& RGPSGravityDriver::PhiFace1(MeshBlock *pmb) const {
  return pmb->ruser_meshblock_data[0];
}

AthenaArray<Real>& RGPSGravityDriver::MgravFace1(MeshBlock *pmb) {
  return pmb->ruser_meshblock_data[1];
}

const AthenaArray<Real>& RGPSGravityDriver::MgravFace1(MeshBlock *pmb) const {
  return pmb->ruser_meshblock_data[1];
}

Real RGPSGravityDriver::CellPhi(MeshBlock *pmb, int i) const {
  const auto &Phi = PhiFace1(pmb);
  const Real rc = pmb->pcoord->x1v(i);
  const Real rm = pmb->pcoord->x1f(i);
  const int ig = GlobalRadialIndex(pmb, i);

  if (ig < 0) {
    const Real r = pmb->pcoord->x1v(i);
    const Real rin = r_face_global_(0);
    const Real bh_mass = GetBlackHoleMass();

    return Phi_face_global_(0)
      + 0.5*std::log((1.0 - 2.0*bh_mass/r)/(1.0 - 2.0*bh_mass/rin));
  }

  if(ig >= 0){
    const Real r = pmb->pcoord->x1v(i);
    const Real mout = mgrav_face_global_(nr_);
    return 0.5*std::log(1.0 - 2.0*mout/r);
  }
  
  const Real dPhi_dr = dPhi_dr_cell_global_(ig);
  return Phi(i) + (rc-rm)*dPhi_dr;
}

Real RGPSGravityDriver::CellMgrav(MeshBlock *pmb, int i) const {
  const int ig = GlobalRadialIndex(pmb, i);
  if (ig < 0) return GetBlackHoleMass();
  if ( ig >= nr_ ){
    return mgrav_face_global_(nr_);
  }
  const Real X_sq  = X_sq_cell_global_(ig);
  const Real r = r_cell_global_(ig);
  const Real mgrav = 0.5*r*(1.0-1.0/X_sq);
  return mgrav;
}

Real RGPSGravityDriver::CellXsq(MeshBlock *pmb, int i) const {
  int ig = GlobalRadialIndex(pmb, i);

  if (ig < 0) {
    const Real mgrav = mgrav_face_global_(0); // = bh_mass
    const Real r = pmb->pcoord->x1v(i);
    return 1.0 / ( 1.0 - 2.0*mgrav/r);
  }
  
  if (ig >= nr_){
    const Real mgrav = mgrav_face_global_(ig+1);
    const Real r = pmb->pcoord->x1v(i);
    return 1.0 / ( 1.0 - 2.0*mgrav/r);
  };
  
  return X_sq_cell_global_(ig);
}

int RGPSGravityDriver::NumModelOutputVariables() const {
  return 3;
}

const char *RGPSGravityDriver::ModelOutputVariableName(int n) const {
  if(n==0)return "Phi";
  if(n==1)return "Mgrav";
  if(n==2)return "q";
  return "";
}

Real RGPSGravityDriver::ModelOutputVariable(
    MeshBlock *pmb,
    int n,
    int k, int j, int i) const{
    if (n == 0) return CellPhi(pmb, i);
    if (n == 1) return CellMgrav(pmb, i);
    if (n == 2) return CellDensitizationFactor(pmb, k,j,i);
    return 0.0;
}

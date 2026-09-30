//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file outflow_boundary_data.cpp
//! \brief Outflow data to inject from the inner boundary.

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../athena.hpp"
#include "../defs.hpp"
#include "outflow_boundary_data.hpp"
#include "hdf5_reader.hpp"
#include "../globals.hpp"

#ifndef HDF5OUTPUT
#error "gr_collapsar with outflow injection requires HDF5 support"
#endif

#ifdef HDF5OUTPUT

#include <hdf5.h>


namespace {
  static constexpr Real c_table_    = 2.99792458e10;  // cm/s
  static constexpr Real G_table_    = 6.6740e-8;     // cgs
  static constexpr Real Msun_table_ = 1.989e33;     // g

  static constexpr Real length_unit_table_ =
    G_table_*Msun_table_/(c_table_*c_table_); // cm

  static constexpr Real time_unit_table_ =
    G_table_ * Msun_table_ /
    (c_table_ * c_table_ * c_table_);

  static constexpr Real rho_unit_table_ = 
    Msun_table_/(length_unit_table_*length_unit_table_*length_unit_table_);

  static constexpr Real press_unit_table_ = 
    rho_unit_table_*c_table_*c_table_;
}

OutflowBoundaryData::OutflowBoundaryData(const std::string &filename){
  
  const char *var_names[NVAR_OUT] = {"rho", "press", "vx", "vy", "vz", "ye", "entropy", "alpha", "h", "hut", "psi", "qe", "rho_star", "ut", "temp", "w"};

  HDF5TableLoader(filename.c_str(), &table_, NVAR_OUT, var_names, "time_lim", "theta_lim");

  table_.GetX1lim(theta_min, theta_max);
  table_.GetX2lim(time_min, time_max);
  
  int nvar_tmp;
  table_.GetSize(nvar_tmp, ntime, ntheta);
  if (Globals::my_rank == 0) {
    std::cout << "Nvar = " << nvar_tmp 
	      << ", Ntime = " << ntime
	      << ", Ntheta = " << ntheta << std::endl;
  }
  
  // Read metadata
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  double tmp;

  hid_t dset = H5Dopen2(file, "r_ext", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
  H5Dclose(dset);
  r_ext_ = static_cast<Real>(tmp);

  dset = H5Dopen2(file, "h_min_global", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
  H5Dclose(dset);
  h_min_global_ = static_cast<Real>(tmp);

  H5Fclose(file);

  if (Globals::my_rank == 0) {
    std::cout << "r_ext = " << r_ext_ << " cm"
	      << ", h_min = " << h_min_global_ << std::endl;
  }
}

OutflowBoundaryData::~OutflowBoundaryData() {
}

OutflowState OutflowBoundaryData::Interpolate(Real time_cgs, Real theta) const {

  Real time_geo = time_cgs/time_unit_table_;
  
  Real q[NVAR_OUT];
  table_.interpolate_all(time_geo, theta, q);

  OutflowState state;
  state.rho     = q[IRHO_OUT]*rho_unit_table_;
  state.press   = q[IPRESS_OUT]*press_unit_table_;
  state.vx      = q[IVX_OUT];
  state.vy      = q[IVY_OUT];
  state.vz      = q[IVZ_OUT];
  state.ye      = q[IYE_OUT];
  state.entropy = q[IENTROPY_OUT];
  
  return state;
}

OutflowState OutflowBoundaryData::GetState(int it, int j) const {
  OutflowState state;

  state.rho          = table_.data(IRHO_OUT,     it, j) *rho_unit_table_;
  state.press        = table_.data(IPRESS_OUT,   it, j) *press_unit_table_;
  state.vx           = table_.data(IVX_OUT,      it, j);
  state.vy           = table_.data(IVY_OUT,      it, j);
  state.vz           = table_.data(IVZ_OUT,      it, j);
  state.ye           = table_.data(IYE_OUT,      it, j);
  state.entropy      = table_.data(IENTROPY_OUT, it, j);
  state.alpha        = table_.data(IALPHA_OUT,   it, j);
  state.h            = table_.data(IENTHALPY_OUT,it, j);
  state.hut          = table_.data(IHUT_OUT,     it, j);
  state.psi          = table_.data(IPSI_OUT,     it, j);
  state.qe           = table_.data(IQE_OUT,      it, j);
  state.rho_star     = table_.data(IRHOSTAR_OUT, it, j) *rho_unit_table_;
  state.temp         = table_.data(ITEMP_OUT,    it, j);
  state.ut           = table_.data(IUT_OUT,      it, j);
  state.w            = table_.data(IW_OUT,       it, j);

  return state;
}

Real OutflowBoundaryData::GetTimeMin() const {
  return time_min*time_unit_table_;
}

Real OutflowBoundaryData::GetTimeMax() const {
  return time_max*time_unit_table_;
}

Real OutflowBoundaryData::GetThetaMin() const {
  return theta_min;
}

Real OutflowBoundaryData::GetThetaMax() const {
  return theta_max;
}

Real OutflowBoundaryData::GetNTime() const {
  return ntime;
}

Real OutflowBoundaryData::GetNTheta() const {
  return ntheta;
}

void OutflowBoundaryData::Analyze() const {
  const Real dtheta = 0.5*M_PI/(Real)ntheta;
  const Real dt = (time_max - time_min)/(Real)ntime;

  Real Mej_geom = 0.0;
  Real Mej_bern = 0.0;
  Real Mej_bind = 0.0;
  Real E_geom = 0.0;
  Real E_bern = 0.0;
  Real E_bind = 0.0;
  for(int it=0; it<ntime; ++it){
    Real t = time_min + dt*(Real)it;

    for(int j=0; j<ntheta; ++j){
      Real theta_c = 0.5*dtheta + dtheta * (Real)j;
      Real theta_d = theta_c - 0.5*dtheta;
      Real theta_u = theta_c + 0.5*dtheta;
      
      OutflowState state = GetState(it, j);
      Real dS = r_ext_*r_ext_ * 2.0*M_PI * (std::cos(theta_d) - std::cos(theta_u))*2.0;
      Real vr = state.vx * std::sin(theta_c) + state.vz * std::cos(theta_c);
      Real dM = dS * state.rho_star * vr*c_table_ * dt*time_unit_table_;

      if(state.ut+1.0 < 0.0 and dM>0.0){
	Mej_geom += dM;
	E_geom += dM*(-state.ut-1.0)*c_table_*c_table_;
      }
      if(state.hut+h_min_global_ < 0.0  and dM>0.0){
	Mej_bern += dM;
	E_bern += dM*(-state.hut-h_min_global_)*c_table_*c_table_;
      }
      Real ebind = state.alpha*state.qe - h_min_global_;
      if(ebind > 0.0 and dM>0.0){
	Mej_bind += dM;
	E_bind += dM*ebind*c_table_*c_table_;
      }
      
      // if(Globals::my_rank==0 and j==ntheta-1){
      // 	std::cout << std::setprecision(5)
      // 		  << t * time_unit_table_<< " "
      // 		  << theta_c << " "
      // 		  << state.rho << " "
      // 		  << state.press << " "
      // 		  << state.vx << " "
      // 		  << state.vy << " "
      // 		  << state.vz << " "
      // 		  << state.ye << " "
      // 		  << state.entropy << " "
      // 		  << state.alpha << " "
      // 		  << state.h << " "
      // 		  << state.hut << " "
      // 		  << state.psi << " "
      // 		  << state.qe << " "
      // 		  << state.rho_star << " "
      // 		  << state.temp << " "
      // 		  << state.ut << " "
      // 		  << state.w << std::endl;
      // }
    }
  }
  std::cout << std::setprecision(5)
	    << " Mej(geom) = " << Mej_geom/Msun_table_ << " Msun"
	    << " Mej(bern) = " << Mej_bern/Msun_table_ << " Msun"
	    << " Mej(bind) = " << Mej_bind/Msun_table_ << " Msun" << std::endl;

  std::cout << std::setprecision(5)
	    << " E(geom) = " << E_geom << " erg"
	    << " E(bern) = " << E_bern << " erg"
	    << " E(bind) = " << E_bind << " erg" << std::endl;

  std::abort();

}

#endif  // HDF5OUTPUT


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
#include <fstream>

#include "../athena.hpp"
#include "../defs.hpp"
#include "outflow_boundary_data.hpp"
#include "hdf5_reader.hpp"
#include "../globals.hpp"

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

  //HDF5TableLoader(filename.c_str(), &table_, NVAR_OUT, var_names, "time_lim", "theta_lim");
  HDF5NonUniformTableLoader(filename.c_str(), &table_, NVAR_OUT, var_names, "time", "theta");

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
  
  // Useful variable to be kept in this class
  time_.NewAthenaArray(ntime);
  theta_c_.NewAthenaArray(ntheta);
  theta_d_.NewAthenaArray(ntheta);
  theta_u_.NewAthenaArray(ntheta);

  dset = H5Dopen2(file, "time", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_.data());
  H5Dclose(dset);

  dset = H5Dopen2(file, "theta", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theta_c_.data());
  H5Dclose(dset);
  
  auto fail = [](const std::string &message) {
    std::stringstream msg;
    msg << "### FATAL ERROR in OutflowBoundaryData\n" << message;
    ATHENA_ERROR(msg);
  };
  
  const htri_t has_d = H5Lexists(file, "theta_d", H5P_DEFAULT);
  const htri_t has_u = H5Lexists(file, "theta_u", H5P_DEFAULT);
  
  if (has_d < 0 || has_u < 0) {
    fail("Could not check angular boundary datasets.");
  }
  if ((has_d > 0) != (has_u > 0)) {
    fail("theta_d and theta_u must either both exist or both be absent.");
  }


  if (has_d > 0) {
    auto read_edge = [&](const char *name, AthenaArray<Real> &array) {
      hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
      if (ds < 0) {
	fail(std::string("Cannot open dataset: ") + name);
      }
      
      hid_t space = H5Dget_space(ds);
      hsize_t dims[1] = {0};
      bool valid = false;
      if (space >= 0 && H5Sget_simple_extent_ndims(space) == 1) {
	valid = H5Sget_simple_extent_dims(space, dims, nullptr) >= 0
	  && dims[0] == static_cast<hsize_t>(ntheta);
      }
      if (space >= 0) H5Sclose(space);
      
      if (!valid) {
	H5Dclose(ds);
	fail(std::string(name) + " must be 1D with length ntheta.");
      }
      
#if SINGLE_PRECISION_ENABLED
      hid_t memory_type = H5T_NATIVE_FLOAT;
#else
      hid_t memory_type = H5T_NATIVE_DOUBLE;
#endif
      herr_t status = H5Dread(ds, memory_type, H5S_ALL, H5S_ALL,
			      H5P_DEFAULT, array.data());
      H5Dclose(ds);
      if (status < 0) {
	fail(std::string("Cannot read dataset: ") + name);
      }
    };
    
    
    read_edge("theta_d", theta_d_);
    read_edge("theta_u", theta_u_);
    
  } else {
    
    // Legacy format: uniform midpoint samples covering [0, pi/2].
    if (ntheta < 2) {
      fail("At least two angular samples are required.");
    }
    
    const Real half_pi = 0.5 * std::acos(-1.0);
    const Real dtheta = half_pi / static_cast<Real>(ntheta);
    
    // Allow rounding in the original ASCII output and in Real.
#if SINGLE_PRECISION_ENABLED
    const Real tolerance = 1.0e-6;
#else
    const Real tolerance = 1.0e-10;
#endif
    
    for (int j = 0; j < ntheta; ++j) {
      const Real expected = (static_cast<Real>(j) + 0.5) * dtheta;
      if (!std::isfinite(theta_c_(j)) || std::abs(theta_c_(j) - expected) > tolerance) {
	fail("Missing angular boundaries: theta is not a uniform "
	     "midpoint grid covering [0, pi/2].");
      }
      theta_d_(j) = static_cast<Real>(j) * dtheta;
      theta_u_(j) = static_cast<Real>(j + 1) * dtheta;
    }
  }
  
  if (Globals::my_rank == 0) {
    std::cout << (has_d > 0
		  ? "Angular boundaries read from HDF5."
		  : "Angular boundaries reconstructed from uniform midpoint samples.")
	      << std::endl;
  }
  // dset = H5Dopen2(file, "theta_d", H5P_DEFAULT);
  // H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theta_d_.data());
  // H5Dclose(dset);

  // dset = H5Dopen2(file, "theta_u", H5P_DEFAULT);
  // H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theta_u_.data());
  // H5Dclose(dset);

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

int OutflowBoundaryData::GetNTime() const {
  return ntime;
}

int OutflowBoundaryData::GetNTheta() const {
  return ntheta;
}

Real OutflowBoundaryData::GetExtractionRadius() const {
  return r_ext_;
}

AveragedOutflowState OutflowBoundaryData::GetAveragedData(int it, Real theta1, Real theta2) const {

  if(it < 0 || ntime <= it){
    std::stringstream msg;
    msg << "### FATAL ERROR in OutflowBoundaryData\n"
	<< "time index out of bound.";
    ATHENA_ERROR(msg);
  }
  
  if( theta1 >= theta2 ){
    std::stringstream msg;
    msg << "### FATAL ERROR in OutflowBoundaryData\n"
	<< "theta1 must be smaller than theta2.";
    ATHENA_ERROR(msg);
  }

  int ilo, ihi;

  ilo=0;
  ihi=ntheta-1;
  while(ihi-ilo>1){
    int i=(ihi+ilo)/2;
    if(theta_d_(i) < theta1){
      ilo=i;
    }else{
      ihi=i;
    }
  }
  // theta_d_(ilo) < theta1 < theta_d_(ihi)
  const int i1 = ilo; 
  
  ilo=0;
  ihi=ntheta-1;
  while(ihi-ilo>1){
    int i=(ihi+ilo)/2;
    if(theta_u_(i) < theta2){
      ilo=i;
    }else{
      ihi=i;
    }
  }
  // theta_u_(ilo) < theta2 < theta_u_(ihi)
  const int i2 = ihi; 

  Real ave_uu1 = 0.0;
  Real ave_uu2 = 0.0;
  Real ave_uu3 = 0.0;
  Real ave_rho = 0.0;
  Real ave_rho_star = 0.0;
  Real ave_ye = 0.0;
  Real ave_entropy = 0.0;
  Real ave_press= 0.0;
  for(int i=i1; i<=i2; ++i){
    Real lo = std::max(theta1, theta_d_(i));
    Real hi = std::min(theta2, theta_u_(i));
    if (hi <= lo) continue;
    const Real dOmega = (std::cos(lo) - std::cos(hi));

    Real costheta = std::cos(theta_c_(i));
    Real sintheta = std::sin(theta_c_(i));

    // Cylindrical (y=0) to orthonormal matrix.
    Real m1_R = sintheta;
    Real m1_z = costheta;

    Real m2_R = costheta;
    Real m2_z = -sintheta;
    
    Real m3_y = 1.0;
    
    const OutflowState state = GetState(it, i);
    
    const Real psi2 = SQR(state.psi);
    // Velocity in orthonormal frame. Assume beta^i = 0.
    Real uu1 = (m1_R*state.vx + m1_z*state.vz) *psi2*state.w/ state.alpha;
    Real uu2 = (m2_R*state.vx + m2_z*state.vz) *psi2*state.w/ state.alpha;
    Real uu3 =  m3_y*state.vy                  *psi2*state.w/ state.alpha;
    Real rho = state.rho;
    Real rho_star = state.rho_star;
    
    ave_rho_star += dOmega*state.rho_star;
    ave_rho += dOmega*rho;
    ave_uu1 += dOmega*rho*uu1;
    ave_uu2 += dOmega*rho*uu2;
    ave_uu3 += dOmega*rho*uu3;
    ave_ye  += dOmega*rho_star*state.ye;
    ave_entropy += dOmega*rho_star*state.entropy;
    ave_press += dOmega*state.press;
  }
  
  Real dOmega_cell = std::cos(theta1)-std::cos(theta2);
  ave_rho_star /= dOmega_cell;
  ave_rho /= dOmega_cell;
  ave_uu1 /= ave_rho*dOmega_cell;
  ave_uu2 /= ave_rho*dOmega_cell;
  ave_uu3 /= ave_rho*dOmega_cell;
  ave_ye  /= ave_rho_star*dOmega_cell;
  ave_entropy /= ave_rho_star*dOmega_cell;
  ave_press /= dOmega_cell;
  
  AveragedOutflowState ave_state;
  ave_state.rho = ave_rho;
  ave_state.press = ave_press;
  ave_state.u1hat = ave_uu1;
  ave_state.u2hat = ave_uu2;
  ave_state.u3hat = ave_uu3;
  ave_state.ye = ave_ye;
  
  return ave_state;
}

AveragedOutflowState OutflowBoundaryData::GetPointData(int it, Real theta) const {

  if(it < 0 || ntime <= it){
    std::stringstream msg;
    msg << "### FATAL ERROR in OutflowBoundaryData\n"
	<< "time index out of bound.";
    ATHENA_ERROR(msg);
  }

  int ilo, ihi;

  ilo=0;
  ihi=ntheta;
  while(ihi-ilo>1){
    int i=(ihi+ilo)/2;
    if(theta_d_(i) <= theta){
      ilo=i;
    }else{
      ihi=i;
    }
  }
  // theta_d_(ilo) < theta < theta_d_(ihi)
  const int i = ilo; 

  const OutflowState state = GetState(it, i);

  Real costheta = std::cos(theta_c_(i));
  Real sintheta = std::sin(theta_c_(i));
  
  // Cylindrical (y=0) to orthonormal matrix.
  Real m1_R = sintheta;
  Real m1_z = costheta;
  
  Real m2_R = costheta;
  Real m2_z = -sintheta;
  
  Real m3_y = 1.0;
    
  const Real psi2 = SQR(state.psi);
  // Velocity in orthonormal frame. Assume beta^i = 0.
  Real uu1 = (m1_R*state.vx + m1_z*state.vz) *psi2*state.w/ state.alpha;
  Real uu2 = (m2_R*state.vx + m2_z*state.vz) *psi2*state.w/ state.alpha;
  Real uu3 =  m3_y*state.vy                  *psi2*state.w/ state.alpha;
  Real rho = state.rho;
  Real rho_star = state.rho_star;
  
  AveragedOutflowState ave_state;
  ave_state.rho = state.rho;
  ave_state.press = state.press;
  ave_state.u1hat = uu1;
  ave_state.u2hat = uu2;
  ave_state.u3hat = uu3;
  ave_state.ye = state.ye;
  
  return ave_state;
}

void OutflowBoundaryData::Analyze(Real gamma_ad) const {
  //const Real dtheta = 0.5*M_PI/(Real)ntheta;
  //const Real dt = (time_max - time_min)/(Real)(ntime-1);
  //const Real dt_cgs = dt*time_unit_table_;
  const Real c2 = c_table_*c_table_;

  std::ofstream ofs("outflow_analysis.dat");
  ofs << "#";
  for (int i=1; i<=25; ++i){
    ofs << std::setw(16) << i;
  }
  ofs << "\n";
  ofs << std::scientific << std::setprecision(8)
      << "#"
      << std::setw(16) << "time"
      << std::setw(16) << "Mdot_geom"
      << std::setw(16) << "Mej_geom"
      << std::setw(16) << "Edot_geom"
      << std::setw(16) << "Eej_geom"
      << std::setw(16) << "Mdot_bern"
      << std::setw(16) << "Mej_bern"
      << std::setw(16) << "Edot_bern"
      << std::setw(16) << "Eej_bern"
      << std::setw(16) << "Mdot_bind"
      << std::setw(16) << "Mej_bind"
      << std::setw(16) << "Edot_bind"
      << std::setw(16) << "Eej_bind"
      << std::setw(16) << "alpha_av"
      << std::setw(16) << "alpha_min"
      << std::setw(16) << "alpha_max"
      << std::setw(16) << "psi_av"
      << std::setw(16) << "psi_min"
      << std::setw(16) << "psi_max"
      << std::setw(16) << "<h(NR)-h_min>"
      << std::setw(16) << "<h-1>"
      << std::setw(16) << "Mdot_in"
      << std::setw(16) << "Mdot_out"
      << std::setw(16) << "Mass_in"
      << std::setw(16) << "Mass_out"
      << "\n";
  
  Real Mass_in = 0.0;
  Real Mass_out= 0.0;

  Real Mej_geom = 0.0;
  Real Mej_bern = 0.0;
  Real Mej_bind = 0.0;
  Real Eej_geom = 0.0;
  Real Eej_bern = 0.0;
  Real Eej_bind = 0.0;
  Real Mdot_in_prev = 0.0, Mdot_out_prev = 0.0;
  Real Mdot_geom_prev = 0.0, Mdot_bern_prev = 0.0, Mdot_bind_prev = 0.0;
  Real Edot_geom_prev = 0.0, Edot_bern_prev = 0.0, Edot_bind_prev = 0.0;
  for(int it=0; it<ntime; ++it){
    // Real t = time_min + dt*(Real)it;
    Real t = time_(it);
    
    Real alpha_av = 0.0;
    Real alpha_min= 1.0e99;
    Real alpha_max= 0.0;
    Real psi_av = 0.0;
    Real psi_min= 1.0e99;
    Real psi_max= 0.0;

    Real Mdot_in = 0.0;
    Real Mdot_out = 0.0;
    Real Mdot_geom = 0.0;
    Real Mdot_bern = 0.0;
    Real Mdot_bind = 0.0;
    Real Edot_geom = 0.0;
    Real Edot_bern = 0.0;
    Real Edot_bind = 0.0;

    Real dhNR_av = 0.0;
    Real dhAth_av = 0.0;
    
    for(int j=0; j<ntheta; ++j){
      // Real theta_c = 0.5*dtheta + dtheta * (Real)j;
      // Real theta_d = theta_c - 0.5*dtheta;
      // Real theta_u = theta_c + 0.5*dtheta;

      Real theta_c = theta_c_(j);
      Real theta_d = theta_d_(j);
      Real theta_u = theta_u_(j);
      
      OutflowState state = GetState(it, j);
      Real dOmega = 2.0*M_PI * (std::cos(theta_d) - std::cos(theta_u))*2.0;
      Real dS = r_ext_*r_ext_ * dOmega;
      Real vr = state.vx * std::sin(theta_c) + state.vz * std::cos(theta_c);
      Real dM_dt_face = dS * state.rho_star * vr*c_table_;
      
      Real egeom = (-state.ut-1.0)*c2;
      Real ebern = (-state.hut-h_min_global_)*c2;
      Real ebind = (state.alpha*state.qe - h_min_global_)*c2;

      if( vr>0.0){
	Mdot_in += dM_dt_face;
      }else{
	Mdot_out-= dM_dt_face;
      }
      
      if( egeom > 0.0 and vr>0.0){
	Mdot_geom += dM_dt_face;
	Edot_geom += dM_dt_face*egeom;
      }
      if( ebern > 0.0 and vr>0.0){
	Mdot_bern += dM_dt_face;
	Edot_bern += dM_dt_face*ebern;
      }

      if( ebind > 0.0 and vr>0.0){
	Mdot_bind += dM_dt_face;
	Edot_bind += dM_dt_face*ebind;

	Real dhNR  = state.h - h_min_global_;
	dhNR_av += dhNR*dM_dt_face;
	
	Real press = state.press;
	Real rho   = state.rho;
	Real dhAth = gamma_ad/(gamma_ad-1.0) * press/rho/(c_table_*c_table_);
	dhAth_av += dhAth*dM_dt_face;
      }
      
      alpha_av += dOmega * state.alpha;
      alpha_min = std::min(alpha_min, state.alpha);
      alpha_max = std::max(alpha_max, state.alpha);

      psi_av += dOmega * state.psi;
      psi_min = std::min(psi_min, state.psi);
      psi_max = std::max(psi_max, state.psi);


    }
    alpha_av /= (4.0*M_PI);
    psi_av /= (4.0*M_PI);

    if (it > 0) {
      const Real dt_cgs = (time_(it) - time_(it-1)) * time_unit_table_;
      Mass_in += 0.5 * (Mdot_in_prev + Mdot_in) * dt_cgs;
      Mass_out += 0.5 * (Mdot_out_prev + Mdot_out) * dt_cgs;

      Mej_geom += 0.5 * (Mdot_geom_prev + Mdot_geom) * dt_cgs;
      Mej_bern += 0.5 * (Mdot_bern_prev + Mdot_bern) * dt_cgs;
      Mej_bind += 0.5 * (Mdot_bind_prev + Mdot_bind) * dt_cgs;
      Eej_geom += 0.5 * (Edot_geom_prev + Edot_geom) * dt_cgs;
      Eej_bern += 0.5 * (Edot_bern_prev + Edot_bern) * dt_cgs;
      Eej_bind += 0.5 * (Edot_bind_prev + Edot_bind) * dt_cgs;
    }
    Mdot_in_prev = Mdot_in;
    Mdot_out_prev = Mdot_out;
    Mdot_geom_prev = Mdot_geom;
    Mdot_bern_prev = Mdot_bern;
    Mdot_bind_prev = Mdot_bind;
    Edot_geom_prev = Edot_geom;
    Edot_bern_prev = Edot_bern;
    Edot_bind_prev = Edot_bind;

    if(Mdot_bind > 0.0){
      dhNR_av /= Mdot_bind;
      dhAth_av /= Mdot_bind;
    }
    
    ofs << std::scientific << std::setprecision(8)
	<< " "
	<< std::setw(16) << t*time_unit_table_
	<< std::setw(16) << Mdot_geom /Msun_table_
	<< std::setw(16) << Mej_geom /Msun_table_
	<< std::setw(16) << Edot_geom 
	<< std::setw(16) << Eej_geom 
	<< std::setw(16) << Mdot_bern /Msun_table_
	<< std::setw(16) << Mej_bern /Msun_table_
	<< std::setw(16) << Edot_bern 
	<< std::setw(16) << Eej_bern 
	<< std::setw(16) << Mdot_bind /Msun_table_
	<< std::setw(16) << Mej_bind /Msun_table_
	<< std::setw(16) << Edot_bind 
	<< std::setw(16) << Eej_bind 
	<< std::setw(16) << alpha_av
	<< std::setw(16) << alpha_min
	<< std::setw(16) << alpha_max
	<< std::setw(16) << psi_av
	<< std::setw(16) << psi_min
	<< std::setw(16) << psi_max
	<< std::setw(16) << dhNR_av
	<< std::setw(16) << dhAth_av
	<< std::setw(16) << Mdot_in /Msun_table_
	<< std::setw(16) << Mdot_out /Msun_table_
	<< std::setw(16) << Mass_in /Msun_table_
	<< std::setw(16) << Mass_out /Msun_table_
	<< "\n";
  }
  std::cout << std::setprecision(5)
	    << " Mej(geom) = " << Mej_geom/Msun_table_ << " Msun"
	    << " Mej(bern) = " << Mej_bern/Msun_table_ << " Msun"
	    << " Mej(bind) = " << Mej_bind/Msun_table_ << " Msun" << std::endl;

  std::cout << std::setprecision(5)
	    << " E(geom) = " << Eej_geom << " erg"
	    << " E(bern) = " << Eej_bern << " erg"
	    << " E(bind) = " << Eej_bind << " erg" << std::endl;

}

#else

namespace {
void HDF5RequiredError() {
  std::stringstream msg;
  msg << "### FATAL ERROR: Outflow injection requires HDF5 support."
      << std::endl;
  ATHENA_ERROR(msg);
}
}


OutflowBoundaryData::OutflowBoundaryData(const std::string &filename) {
  std::stringstream msg;
  msg << "### FATAL ERROR in OutflowBoundaryData" << std::endl
      << "Outflow injection requires HDF5 support." << std::endl
      << "Reconfigure Athena++ with -hdf5." << std::endl;
  ATHENA_ERROR(msg);
}

OutflowBoundaryData::~OutflowBoundaryData() {
}

OutflowState OutflowBoundaryData::Interpolate(
    Real time, Real theta) const {
  HDF5RequiredError();
  return OutflowState{};
}

OutflowState OutflowBoundaryData::GetState(int it, int j) const {
  HDF5RequiredError();
  return OutflowState{};
}

Real OutflowBoundaryData::GetTimeMin() const {
  HDF5RequiredError();
  return 0.0;
}

Real OutflowBoundaryData::GetTimeMax() const {
  HDF5RequiredError();
  return 0.0;
}

Real OutflowBoundaryData::GetThetaMin() const {
  HDF5RequiredError();
  return 0.0;
}

Real OutflowBoundaryData::GetThetaMax() const {
  HDF5RequiredError();
  return 0.0;
}

int OutflowBoundaryData::GetNTime() const {
  HDF5RequiredError();
  return 0;
}

int OutflowBoundaryData::GetNTheta() const {
  HDF5RequiredError();
  return 0;
}

void OutflowBoundaryData::Analyze(Real gamma_ad) const {
  HDF5RequiredError();
  return;
}

Real OutflowBoundaryData::GetExtractionRadius() const {
  HDF5RequiredError();
  return 0.0;
}

#endif  // HDF5OUTPUT


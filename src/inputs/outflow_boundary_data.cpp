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

#include "../athena.hpp"
#include "../defs.hpp"
#include "outflow_boundary_data.hpp"
#include "hdf5_reader.hpp"

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
  
  const char *var_names[NVAR_OUT] = {"rho", "press", "vx", "vy", "vz", "ye", "entropy"};

  HDF5TableLoader(filename.c_str(), &table_, NVAR_OUT, var_names, "time_lim", "theta_lim");

  table_.GetX1lim(theta_min, theta_max);
  table_.GetX2lim(time_min, time_max);
  
  int nvar_tmp;
  table_.GetSize(nvar_tmp, ntime, ntheta);
  std::cout << "Nvar = " << nvar_tmp 
    << ", Ntime = " << ntime
    << ", Ntheta = " << ntheta << std::endl;
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

#endif  // HDF5OUTPUT


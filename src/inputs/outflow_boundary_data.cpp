//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file progenitor_reader.cpp
//! \brief Reader for standardized 1D progenitor profiles stored in HDF5.

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "../athena.hpp"
#include "../defs.hpp"
#include "outflow_boundary_data.hpp"
#include "hdf5_reader.hpp"

#ifndef HDF5OUTPUT
#error "gr_collapsar with outflow injection requires HDF5 support"
#endif

#ifdef HDF5OUTPUT

#include <hdf5.h>

OutflowBoundaryData::OutflowBoundaryData(const std::string &filename){
  
  const char *var_names[NVAR_OUT] = {"rho", "press", "vx", "vy", "vz", "ye", "entropy"};

  HDF5TableLoader(filename.c_str(), &table_, NVAR_OUT, var_names, "time_lim", "theta_lim");
}

OutflowBoundaryData::~OutflowBoundaryData() {
}

OutflowState OutflowBoundaryData::Interpolate(Real time, Real theta) const {

  Real q[NVAR_OUT];
  table_.interpolate_all(time, theta, q);

  OutflowState state;
  state.rho     = q[IRHO_OUT];
  state.press   = q[IPRESS_OUT];
  state.vx      = q[IVX_OUT];
  state.vy      = q[IVY_OUT];
  state.vz      = q[IVZ_OUT];
  state.ye      = q[IYE_OUT];
  state.entropy = q[IENTROPY_OUT];
  
  return state;
}

#endif  // HDF5OUTPUT

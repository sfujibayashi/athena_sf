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

#ifdef HDF5OUTPUT

#include <hdf5.h>

OutflowBoundaryData::OutflowBoundaryData(const std::string &filename){
}

OutflowBoundaryData::~OutflowBoundaryData() {
}

Real OutflowBoundaryData::InterpolateRho(Real time, Real theta) const {
  return 0.0;
}

Real OutflowBoundaryData::InterpolatePress(Real time, Real theta) const {
  return 0.0;
}

#endif  // HDF5OUTPUT

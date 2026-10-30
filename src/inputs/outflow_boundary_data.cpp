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
#include <algorithm>
#include <limits>

#include "../athena.hpp"
#include "../defs.hpp"
#include "outflow_boundary_data.hpp"
#include "hdf5_reader.hpp"
#include "../globals.hpp"
#include "../utils/utils.hpp"

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

void GaussianError(const std::string &message) {
  std::stringstream msg;
  msg << "### FATAL ERROR in outflow Gaussian smoothing\n" << message;
  ATHENA_ERROR(msg);
}

void GaussianLegendreRule(int n, std::vector<Real> &x,
                          std::vector<Real> &w) {
  if (n < 4 || n > 32) GaussianError("Quadrature order must be between 4 and 32.");
  AthenaArray<Real> nodes, weights;
  nodes.NewAthenaArray(n);
  weights.NewAthenaArray(n);
  Gauleg(n, -1.0, 1.0, nodes, weights);
  x.assign(nodes.data(), nodes.data()+n);
  w.assign(weights.data(), weights.data()+n);
}

// P_l(cos(theta)) and B_l = sin(theta) P_l'(cos(theta)). The associated
// Legendre recurrence avoids dividing by sin(theta) near the polar axis.
void GaussianBasis(Real theta, int lmax, std::vector<Real> &p,
                    std::vector<Real> &b) {
  const Real half_pi = 0.5*std::acos(-1.0);
  const Real mu = theta == half_pi ? 0.0 : std::cos(theta);
  const Real st = std::sin(theta);
  p.resize(lmax+1);
  b.resize(lmax+1);
  p[0] = 1.0;
  b[0] = 0.0;
  p[1] = mu;
  b[1] = st;
  for (int l = 2; l <= lmax; ++l) {
    p[l] = ((2*l-1)*mu*p[l-1]-(l-1)*p[l-2])/l;
    b[l] = ((2*l-1)*mu*b[l-1]-l*b[l-2])/(l-1);
  }
}

struct GaussianAngularWeights {
  std::vector<Real> scalar, polar, azimuthal;
};

GaussianAngularWeights MakeGaussianWeights(
    const AthenaArray<Real> &lower, const AthenaArray<Real> &upper,
    int nsource, const AthenaArray<Real> &centers, int js, int je,
    Real fwhm_deg, int order) {

  const Real pi = std::acos(-1.0), half_pi = pi/2.0;

  if (!std::isfinite(fwhm_deg) || fwhm_deg <= 0.0 || fwhm_deg > 180.0)
    GaussianError("FWHM must be finite and in (0, 180] degrees.");
  if (js < 0 || je < js || je >= centers.GetSize())
    GaussianError("Target angular index is out of bounds.");

  const Real sigma = fwhm_deg*pi/180.0/std::sqrt(8.0*std::log(2.0));
  // Bound the discarded heat-kernel tail, including its mode multiplicity.
  const Real required = std::sqrt(-2.0*std::log(1.0e-14*sigma*sigma))/sigma;
  if (!std::isfinite(required) || required > 8190.0)
    GaussianError("FWHM is too small for the spectral cutoff (8192 modes).");

  const int lmax = static_cast<int>(std::ceil(required))+2;
  const int ns = lmax/2+1, nj = je-js+1;
  const Real angle_tol = 128.0*std::numeric_limits<Real>::epsilon()*half_pi;
  std::vector<Real> nodes, quad;

  GaussianLegendreRule(order, nodes, quad);
  std::vector<Real> scalar_integral(static_cast<std::size_t>(ns)*nsource, 0.0);
  std::vector<Real> vector_integral(static_cast<std::size_t>(lmax+1)*nsource, 0.0);
  std::vector<Real> p, b, si(lmax+1), vi(lmax+1);
  if (nsource < 2 || std::abs(lower(0)) > angle_tol
      || std::abs(upper(nsource-1)-half_pi) > angle_tol)
    GaussianError("Source angular intervals must cover [0, pi/2].");
  for (int i = 0; i < nsource; ++i) {
    if (!std::isfinite(lower(i)) || !std::isfinite(upper(i))
        || upper(i) <= lower(i)
        || (i > 0 && std::abs(lower(i)-upper(i-1)) > angle_tol))
      GaussianError("Source angular intervals must be contiguous and ordered.");
    const Real lo = std::max(0.0, static_cast<Real>(lower(i)));
    const Real hi = std::min(half_pi, static_cast<Real>(upper(i)));
    std::fill(si.begin(), si.end(), 0.0);
    std::fill(vi.begin(), vi.end(), 0.0);
    // Subdivide unusually wide source intervals, not the target Athena grid.
    const int pieces = std::max(1, static_cast<int>(std::ceil((hi-lo)/sigma)));
    for (int piece = 0; piece < pieces; ++piece) {
      const Real a = lo+(hi-lo)*piece/pieces;
      const Real c = lo+(hi-lo)*(piece+1)/pieces;
      const Real middle = 0.5*(a+c), half = 0.5*(c-a);
      for (int k = 0; k < order; ++k) {
        const Real theta = middle+half*nodes[k];
        GaussianBasis(theta, lmax, p, b);
        const Real weight = half*quad[k]*std::sin(theta);
        for (int l = 0; l <= lmax; l += 2) si[l] += weight*p[l];
        for (int l = 1; l <= lmax; ++l) vi[l] += weight*b[l];
      }
    }
    for (int l = 0; l <= lmax; ++l) {
      const Real damping = std::exp(-0.5*l*(l+1)*sigma*sigma);
      if (l%2 == 0)
        scalar_integral[static_cast<std::size_t>(l/2)*nsource+i]
            = (2*l+1)*damping*si[l];
      if (l > 0)
        vector_integral[static_cast<std::size_t>(l)*nsource+i]
            = (2*l+1)*damping*vi[l]/(l*(l+1));
    }
  }
  GaussianAngularWeights result;
  const std::size_t size = static_cast<std::size_t>(nj)*nsource;
  result.scalar.assign(size, 0.0);
  result.polar.assign(size, 0.0);
  result.azimuthal.assign(size, 0.0);
  for (int j = js; j <= je; ++j) {
    const Real theta = centers(j);
    if (!std::isfinite(theta) || theta < -angle_tol || theta > half_pi+angle_tol)
      GaussianError("Target angle must lie in [0, pi/2].");
    GaussianBasis(std::max(0.0, std::min(half_pi, theta)), lmax, p, b);
    const std::size_t offset = static_cast<std::size_t>(j-js)*nsource;
    for (int l = 0; l <= lmax; ++l) {
      if (l%2 == 0) {
        const std::size_t row = static_cast<std::size_t>(l/2)*nsource;
        for (int i = 0; i < nsource; ++i)
          result.scalar[offset+i] += p[l]*scalar_integral[row+i];
      }
      if (l > 0) {
        const std::size_t row = static_cast<std::size_t>(l)*nsource;
        std::vector<Real> &weights = l%2 == 0 ? result.polar : result.azimuthal;
        for (int i = 0; i < nsource; ++i)
          weights[offset+i] += b[l]*vector_integral[row+i];
      }
    }
    Real sum = 0.0;
    for (int i = 0; i < nsource; ++i) {
      const std::size_t index = offset+i;
      if (result.scalar[index] < -1.0e-10 || result.polar[index] < -1.0e-10
          || result.azimuthal[index] < -1.0e-10)
        GaussianError("Kernel quadrature is inaccurate; increase quadrature order.");
      // The exact heat kernels are nonnegative. Remove spectral roundoff tails.
      result.scalar[index] = std::max(0.0, result.scalar[index]);
      result.polar[index] = std::max(0.0, result.polar[index]);
      result.azimuthal[index] = std::max(0.0, result.azimuthal[index]);
      sum += result.scalar[index];
    }
    if (!std::isfinite(sum) || std::abs(sum-1.0) > 1.0e-8)
      GaussianError("Scalar kernel normalization failed; increase quadrature order.");
    for (int i = 0; i < nsource; ++i) result.scalar[offset+i] /= sum;
    // Do not normalize vector kernels: axis/equator parity attenuates them.
  }
  return result;
}
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

  Real tmp;

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

Real OutflowBoundaryData::GetHMin() const {
  return h_min_global_;
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
  Real ave_press= 0.0; // <P>
  Real ave_E = 0.0;  // <rho*h*W^2 - P>
  Real ave_alpha = 0.0;
  Real ave_psi = 0.0;
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

    ave_E   += dOmega*state.qe*state.rho*state.w;
    ave_alpha += dOmega*state.alpha;
    ave_psi += dOmega*state.psi;
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
  ave_E     /= dOmega_cell;
  ave_alpha /= dOmega_cell;
  ave_psi   /= dOmega_cell;
  
  AveragedOutflowState ave_state;
  // primitive
  ave_state.rho = ave_rho;
  ave_state.press = ave_press;
  ave_state.u1hat = ave_uu1;
  ave_state.u2hat = ave_uu2;
  ave_state.u3hat = ave_uu3;
  ave_state.ye = ave_ye;
  ave_state.E = ave_E*SQR(c_table_);
  ave_state.alpha = ave_alpha;
  ave_state.psi = ave_psi;
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
  ave_state.E = state.qe*state.rho*state.w*SQR(c_table_);
  ave_state.alpha = state.alpha;
  ave_state.psi = state.psi;
  
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

Real OutflowBoundaryData::GetTimeCGS(int it) const {
  if(it<0 || ntime<=it){
    std::stringstream msg;
    msg << "### FATAL ERROR in OutflowBoundaryData\n"
	<< "it out of bound: it=" << it
	<< "ntime = " << ntime;
    ATHENA_ERROR(msg);
  }
  return time_(it)*time_unit_table_;
}


// OutflowBoundaryTable
void OutflowBoundaryData::FillGaussianTable(
    const AthenaArray<Real> &theta_centers, int js, int je,
    Real fwhm_deg, int quadrature_order, AthenaArray<Real> &out) const {
  const auto weights = MakeGaussianWeights(theta_d_, theta_u_, ntheta,
      theta_centers, js, je, fwhm_deg, quadrature_order);
  
  // Source values are assumed to be piecewise constant in their original angular intervals.
  // Transform at the original sample angle, exactly as in GetAveragedData.
  std::vector<Real> q(static_cast<std::size_t>(ntheta)*10);
  for (int it = 0; it < ntime; ++it) {
    for (int i = 0; i < ntheta; ++i) {
      const OutflowState s = GetState(it, i);
      const Real st = std::sin(theta_c_(i)), ct = std::cos(theta_c_(i));
      const Real factor = s.psi*s.psi*s.w/s.alpha;
      Real *v = &q[static_cast<std::size_t>(i)*10];
      v[0] = s.rho;
      v[1] = s.press;
      v[2] = s.rho*(s.vx*st+s.vz*ct)*factor;
      v[3] = s.rho*(s.vx*ct-s.vz*st)*factor;
      v[4] = s.rho*s.vy*factor;
      v[5] = s.rho_star*s.ye;
      v[6] = s.rho_star;
      v[7] = s.qe*s.rho*s.w*c_table_*c_table_;
      v[8] = s.alpha;
      v[9] = s.psi;
    }
    for (int j = js; j <= je; ++j) {
      Real v[10] = {};
      const std::size_t offset = static_cast<std::size_t>(j-js)*ntheta;
      for (int i = 0; i < ntheta; ++i) {
        const Real ws = weights.scalar[offset+i];
        const Real wt = weights.polar[offset+i];
        const Real wp = weights.azimuthal[offset+i];
        const Real *s = &q[static_cast<std::size_t>(i)*10];
        v[0] += ws*s[0];
        v[1] += ws*s[1];
        v[2] += ws*s[2];
        v[3] += wt*s[3];
        v[4] += wp*s[4];
        for (int k = 5; k < 10; ++k) v[k] += ws*s[k];
      }
      if (!(v[0] > 0.0) || !(v[6] > 0.0))
        GaussianError("Smoothed rho and rho_star must be positive.");
      out(it,j,IRHO_AVE) = v[0];
      out(it,j,IPRESS_AVE) = v[1];
      out(it,j,IU1HAT_AVE) = v[2]/v[0];
      out(it,j,IU2HAT_AVE) = v[3]/v[0];
      out(it,j,IU3HAT_AVE) = v[4]/v[0];
      out(it,j,IYE_AVE) = v[5]/v[6];
      out(it,j,IE_AVE) = v[7];
      out(it,j,IALPHA_AVE) = v[8];
      out(it,j,IPSI_AVE) = v[9];
    }
  }
}

// constructor
OutflowBoundaryTable::OutflowBoundaryTable(const OutflowBoundaryData &source,
					   const AthenaArray<Real> &theta_faces,
					   const int js, const int je)
    : OutflowBoundaryTable(source, theta_faces, theta_faces, js, je, 0.0, 8) {
}

OutflowBoundaryTable::OutflowBoundaryTable(const OutflowBoundaryData &source,
    const AthenaArray<Real> &theta_faces,
    const AthenaArray<Real> &theta_centers, int js, int je,
    Real smoothing_fwhm_deg, int smoothing_quadrature_order) {
  if (!std::isfinite(smoothing_fwhm_deg) || smoothing_fwhm_deg < 0.0)
    GaussianError("Smoothing FWHM must be finite and nonnegative.");
  //
  r_ext_ = source.GetExtractionRadius();
  h_min_global_ = source.GetHMin();

  // allocate data
  const int ntheta = theta_faces.GetSize() - 1;
  const int ntime = source.GetNTime();
  ntime_=ntime;
  data_.NewAthenaArray(ntime,ntheta,NVAR_AVE);
  time_.NewAthenaArray(ntime_);

  // input time
  for(int it=0; it<ntime_; ++it){
    time_(it) = source.GetTimeCGS(it);
  }

  if (smoothing_fwhm_deg > 0.0) {
    source.FillGaussianTable(theta_centers, js, je, smoothing_fwhm_deg,
                             smoothing_quadrature_order, data_);
    return;
  }
  
  // average between theta faces
  for(int it=0; it<ntime_; ++it){
    for(int j=js; j<=je; ++j){
      AveragedOutflowState ave_state = source.GetAveragedData(it, theta_faces(j), theta_faces(j+1));
      data_(it,j,IRHO_AVE) = ave_state.rho;
      data_(it,j,IPRESS_AVE) = ave_state.press;
      data_(it,j,IU1HAT_AVE) = ave_state.u1hat;
      data_(it,j,IU2HAT_AVE) = ave_state.u2hat;
      data_(it,j,IU3HAT_AVE) = ave_state.u3hat;
      data_(it,j,IYE_AVE) = ave_state.ye;
      data_(it,j,IE_AVE) = ave_state.E;
      data_(it,j,IALPHA_AVE) = ave_state.alpha;
      data_(it,j,IPSI_AVE) = ave_state.psi;
    }
  }
}

void OutflowBoundaryTable::Analyze(Real gamma_ad,
				   const AthenaArray<Real> &theta,
				   const AthenaArray<Real> &theta_faces,
				   const int js, const int je,
				   AthenaArray<Real> &out) const {

  const Real c2 = c_table_*c_table_;

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
  for(int it=0; it<ntime_; ++it){

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
    
    for(int j=js; j<=je; ++j){

      Real theta_c = theta(j);
      Real theta_d = theta_faces(j);
      Real theta_u = theta_faces(j+1);

      const Real rho   = data_(it,j,IRHO_AVE);
      const Real press = data_(it,j,IPRESS_AVE);      
      const Real u1hat = data_(it,j,IU1HAT_AVE);
      const Real u2hat = data_(it,j,IU2HAT_AVE);
      const Real u3hat = data_(it,j,IU3HAT_AVE);
      const Real E     = data_(it,j,IE_AVE);
      const Real alpha = data_(it,j,IALPHA_AVE);
      const Real psi   = data_(it,j,IPSI_AVE);

      const Real psi2 = psi*psi;
      const Real psi6 = psi2*psi2*psi2;
      const Real ur = u1hat/psi2;
      const Real W = std::sqrt(1.0 + SQR(u1hat) + SQR(u2hat) + SQR(u3hat));
      const Real qe = E/(rho*W*c2);
      
      Real dOmega = 2.0*M_PI * (std::cos(theta_d) - std::cos(theta_u))*2.0;
      Real dS = r_ext_*r_ext_ * dOmega;
      Real dM_dt_face = dS * psi6*alpha * rho * ur*c_table_;
      Real ut = -alpha*W;
      Real h = (E+press)/(rho*SQR(W)*c2);
      
      Real egeom = (-ut-1.0)*c2;
      Real ebern = (-h*ut-h_min_global_)*c2;
      Real ebind = (alpha*qe - h_min_global_)*c2;

      if( ur>0.0){
	Mdot_in += dM_dt_face;
      }else{
	Mdot_out-= dM_dt_face;
      }
      
      if( egeom > 0.0 and ur>0.0){
	Mdot_geom += dM_dt_face;
	Edot_geom += dM_dt_face*egeom;
      }
      if( ebern > 0.0 and ur>0.0){
	Mdot_bern += dM_dt_face;
	Edot_bern += dM_dt_face*ebern;
      }

      if( ebind > 0.0 and ur>0.0){
	Mdot_bind += dM_dt_face;
	Edot_bind += dM_dt_face*ebind;

	Real dhNR  = h - h_min_global_;
	dhNR_av += dhNR*dM_dt_face;
	
	Real dhAth = gamma_ad/(gamma_ad-1.0) * press/(rho*c2);
	dhAth_av += dhAth*dM_dt_face;
      }
      
      alpha_av += dOmega * alpha;
      alpha_min = std::min(alpha_min, alpha);
      alpha_max = std::max(alpha_max, alpha);

      psi_av += dOmega * psi;
      psi_min = std::min(psi_min, psi);
      psi_max = std::max(psi_max, psi);


    }
    alpha_av /= (4.0*M_PI);
    psi_av /= (4.0*M_PI);

    if (it > 0) {
      const Real dt_cgs = time_(it) - time_(it-1);
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

    // store in array
    out(it,IMASS_IN)  += Mass_in  /Msun_table_;
    out(it,IMASS_OUT) += Mass_out /Msun_table_;
    out(it,IMEJ_GEOM) += Mej_geom /Msun_table_;
    out(it,IMEJ_BERN) += Mej_bern /Msun_table_;
    out(it,IMEJ_BIND) += Mej_bind /Msun_table_;
    out(it,IEEJ_GEOM) += Eej_geom;
    out(it,IEEJ_BERN) += Eej_bern;
    out(it,IEEJ_BIND) += Eej_bind;
    out(it,IALPHA_AV) += alpha_av*4.0*M_PI;
    out(it,IALPHA_MIN) = std::min(out(it,IALPHA_MIN), alpha_min);
    out(it,IALPHA_MAX) = std::max(out(it,IALPHA_MAX), alpha_max);
    out(it,IPSI_AV) += psi_av*4.0*M_PI;
    out(it,IPSI_MIN) = std::min(out(it,IPSI_MIN), psi_min);
    out(it,IPSI_MAX) = std::max(out(it,IPSI_MAX), psi_max);
    out(it,IMDOT_IN) += Mdot_in /Msun_table_;
    out(it,IMDOT_OUT) += Mdot_out /Msun_table_;
    out(it,IMDOT_GEOM) += Mdot_geom /Msun_table_;
    out(it,IMDOT_BERN) += Mdot_bern /Msun_table_;
    out(it,IMDOT_BIND) += Mdot_bind /Msun_table_;
    out(it,IEDOT_GEOM) += Edot_geom;
    out(it,IEDOT_BERN) += Edot_bern;
    out(it,IEDOT_BIND) += Edot_bind;
    out(it,IHNR_AV)  += dhNR_av * (Mdot_bind/Msun_table_);
    out(it,IHATH_AV) += dhAth_av * (Mdot_bind/Msun_table_);

  }
}

int OutflowBoundaryTable::GetNTime() const {
  return ntime_;
}

Real OutflowBoundaryTable::GetTime(const int it) const {
  return time_(it);
}

Real OutflowBoundaryTable::GetTimeMax() const {
  return time_(ntime_-1);
}

AveragedOutflowState OutflowBoundaryTable::InterpolatePoint(Real time_cgs, int j) const {
  AveragedOutflowState ave_state;
  
  // find time index
  if(time_cgs <= time_(0)){
    ave_state.rho = data_(0,j,IRHO_AVE);
    ave_state.press = data_(0,j,IPRESS_AVE);
    ave_state.u1hat = data_(0,j,IU1HAT_AVE);
    ave_state.u2hat = data_(0,j,IU2HAT_AVE);
    ave_state.u3hat = data_(0,j,IU3HAT_AVE);
    ave_state.ye = data_(0,j,IYE_AVE);
    return ave_state;
  }

  if(time_cgs >= time_(ntime_-1)){
    ave_state.rho = data_(ntime_-1,j,IRHO_AVE);
    ave_state.press = data_(ntime_-1,j,IPRESS_AVE);
    ave_state.u1hat = data_(ntime_-1,j,IU1HAT_AVE);
    ave_state.u2hat = data_(ntime_-1,j,IU2HAT_AVE);
    ave_state.u3hat = data_(ntime_-1,j,IU3HAT_AVE);
    ave_state.ye = data_(ntime_-1,j,IYE_AVE);
    return ave_state;
  }

  int ilo=0, ihi=ntime_-1;
  while(ihi-ilo>1){
    int i=(ihi+ilo)/2;
    if(time_(i) <= time_cgs){
      ilo=i;
    }else{
      ihi=i;
    }
  }
  const Real w_time = (time_(ihi)-time_cgs)/(time_(ihi)-time_(ilo));
  // interpolate rho*u and rho
  ave_state.rho = w_time*data_(ilo,j,IRHO_AVE) + (1.0-w_time)*data_(ihi,j,IRHO_AVE);
  ave_state.press = w_time*data_(ilo,j,IPRESS_AVE) + (1.0-w_time)*data_(ihi,j,IPRESS_AVE);
  ave_state.u1hat = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU1HAT_AVE) 
		     + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU1HAT_AVE))/ave_state.rho;
  ave_state.u2hat = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU2HAT_AVE) 
		     + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU2HAT_AVE))/ave_state.rho;
  ave_state.u3hat = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU3HAT_AVE) 
		     + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU3HAT_AVE))/ave_state.rho;
  ave_state.ye = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IYE_AVE) 
		  + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IYE_AVE))/ave_state.rho;
  
  return ave_state;
}

void OutflowBoundaryTable::Interpolate(Real time_cgs, int js, int je, AthenaArray<Real> &out) const {
  
  // find time index
  if(time_cgs <= time_(0)){
    for(int j=js; j<=je; ++j){
      for(int n=0; n<NVAR_AVE; ++n){
	out(j,n) = data_(0,j,n);
      }
    }
    return;
  }

  if(time_cgs >= time_(ntime_-1)){
    for(int j=js; j<=je; ++j){
      for(int n=0; n<NVAR_AVE; ++n){
	out(j,n) = data_(ntime_-1,j,n);
      }
    }
    return;
  }

  int ilo=0, ihi=ntime_-1;
  while(ihi-ilo>1){
    int i=(ihi+ilo)/2;
    if(time_(i) <= time_cgs){
      ilo=i;
    }else{
      ihi=i;
    }
  }

  const Real w_time = (time_(ihi)-time_cgs)/(time_(ihi)-time_(ilo));
  for(int j=js; j<=je; ++j){
    Real rho = w_time*data_(ilo,j,IRHO_AVE) + (1.0-w_time)*data_(ihi,j,IRHO_AVE);
    out(j,IRHO_AVE)   = rho;
    out(j,IPRESS_AVE) = w_time*data_(ilo,j,IPRESS_AVE) + (1.0-w_time)*data_(ihi,j,IPRESS_AVE);
    // interpolate rho*u and rho for velocity
    out(j,IU1HAT_AVE) = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU1HAT_AVE) 
			 + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU1HAT_AVE))/rho;
    out(j,IU2HAT_AVE) = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU2HAT_AVE) 
			 + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU2HAT_AVE))/rho;
    out(j,IU3HAT_AVE) = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IU3HAT_AVE) 
			 + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IU3HAT_AVE))/rho;
    out(j,IYE_AVE)    = (w_time*data_(ilo,j,IRHO_AVE)*data_(ilo,j,IYE_AVE) 
			 + (1.0-w_time)*data_(ihi,j,IRHO_AVE)*data_(ihi,j,IYE_AVE))/rho;
  }
  return;
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


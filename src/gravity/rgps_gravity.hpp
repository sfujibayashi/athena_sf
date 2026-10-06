#ifndef GRAVITY_RGPS_GRAVITY_HPP_
#define GRAVITY_RGPS_GRAVITY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monopole_gravity.hpp
//! \brief defines MonopoleGravityDriver class

// C headers

// C++ headers

//
#include "../athena.hpp"

//! \class MonopoleGravity
//! \brief Constructs global radial monopole gravity profiles
class RGPSGravity {
public:
  RGPSGravity() = default;
  ~RGPSGravity() = default;
  
  void ConstructCovariantMetric(Real r, Real theta, Real phi, Real Phi, Real mgrav,
				Real &g00, Real &g01, Real &g02, Real &g03,
				Real &g11, Real &g12, Real &g13,
				Real &g22, Real &g23, Real &g33) const;

  void MetricRadialDerivatives(Real r, Real theta, Real phi, Real Phi, Real mgrav,
			       Real dPhi_dr, Real dmg_dr,
			       Real &d1_g00, Real &d1_g11) const;

};
#endif // GRAVITY_RGPS_GRAVITY_HPP_

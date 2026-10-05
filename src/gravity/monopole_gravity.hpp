#ifndef GRAVITY_MONOPOLE_GRAVITY_HPP_
#define GRAVITY_MONOPOLE_GRAVITY_HPP_
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
class MonopoleGravity {
public:
  MonopoleGravity();
  ~MonopoleGravity();
  
  void ConstructCovariantMetric(Real r, Real theta, Real phi, Real Psi, Real dm, Real bh_mass,
				Real &g00, Real &g01, Real &g02, Real &g03,
				Real &g11, Real &g12, Real &g13,
				Real &g22, Real &g23, Real &g33) const;
};

#endif // GRAVITY_MONOPOLE_GRAVITY_HPP_

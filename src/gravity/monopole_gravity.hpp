#ifndef GRAVITY_MONOPOLE_GRAVITY_HPP_
#define GRAVITY_MONOPOLE_GRAVITY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monopole_gravity.hpp
//! \brief defines MonopoleGravity class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class Mesh;
class ParameterInput;

//! \class MonopoleGravity
//! \brief Constructs global radial monopole gravity profiles
class MonopoleGravity {
public:
  MonopoleGravity(Mesh *pm, ParameterInput *pin);
  ~MonopoleGravity();
  
  void Update();
  
private:
  Mesh *pmy_mesh_;

  int nr_;

  AthenaArray<Real> dm_shell_global_;
  AthenaArray<Real> delta_m_face_global_;
  AthenaArray<Real> psi_face_global_;

};

#endif // GRAVITY_MONOPOLE_GRAVITY_HPP_

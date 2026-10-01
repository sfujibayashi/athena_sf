#ifndef INPUTS_PROGENITOR_READER_HPP_
#define INPUTS_PROGENITOR_READER_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file progenitor_reader.hpp
//! \brief Reader for standardized 1D progenitor profiles stored in HDF5.
//!
//! Expected datasets:
//!
//!
//!   /mass_face    [N+1]  enclosed baryon mass [g]
//!   /radius_face  [N+1]  radius [cm]
//!
//!   /mass       [N]  enclosed baryon mass [g]
//!   /radius     [N]  radius [cm]
//!   /rho        [N]  rest-mass density [g cm^-3]
//!   /press      [N]  pressure [dyn cm^-2]
//!   /csound     [N]  sound speed [cm s^-1]
//!
//! At least one of:
//!
//!   /jrot       [N]  specific angular momentum [cm^2 s^-1]
//!   /omega      [N]  angular velocity [s^-1]
//!
//! Optional:
//!
//!   /ye         [N]
//!   /temp       [N]  temperature [K]
//!
//! All arrays must use center-to-surface ordering and have the same length.

#include <cstddef>
#include <string>
#include <vector>
#include <array>

#include "../athena.hpp"

enum ProgenitorSpecies {
  IPROG_NEUT = 0,
  IPROG_PROT,
  IPROG_H1,
  IPROG_HE3,
  IPROG_HE4,
  IPROG_C12,
  IPROG_N14,
  IPROG_O16,
  IPROG_NE20,
  IPROG_MG24,
  IPROG_SI28,
  IPROG_S32,
  IPROG_AR36,
  IPROG_CA40,
  IPROG_TI44,
  IPROG_CR48,
  IPROG_CR56,
  IPROG_FE52,
  IPROG_FE54,
  IPROG_FE56,
  IPROG_NI56,
  NPROG_SPECIES
};

struct ProgenitorProfile {
  std::vector<Real> mass_face;
  std::vector<Real> radius_face;

  std::vector<Real> mass;
  std::vector<Real> radius;
  std::vector<Real> rho;
  std::vector<Real> press;
  std::vector<Real> csound;

  std::vector<Real> jrot;
  std::vector<Real> omega;

  std::vector<Real> ye;
  std::vector<Real> temp;

  std::array<std::vector<Real>, NPROG_SPECIES> x;

  bool has_jrot = false;
  bool has_omega = false;
  bool has_ye = false;
  bool has_temp = false;
  bool has_composition = false;

  int ncell;
  int nface; // ncell + 1

  std::size_t size() const {
    return ncell;
  }
};

ProgenitorProfile ReadProgenitorProfile(const std::string &filename);

#endif  // INPUTS_PROGENITOR_READER_HPP_

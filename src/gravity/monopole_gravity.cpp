
// C++ headers
#include <cmath>
#include <iostream>

//
#include "monopole_gravity.hpp"
#include "../athena.hpp"

void MonopoleGravity::ConstructCovariantMetric(
    Real r, Real theta, Real phi, Real Psi, Real dm, Real bh_mass,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const {

  const Real sintheta = std::sin(theta);
  const Real f = 1.0 - 2.0*bh_mass/r;

  g00 = -f + 2.0*dm/r + 2.0*f*Psi;
  g01 = 0.0;
  g02 = 0.0;
  g03 = 0.0;
  
  g11 = 1.0/f + 2.0*dm/(r*f*f);
  g12 = 0.0;
  g13 = 0.0;
  
  g22 = r*r;
  g23 = 0.0;
  g33 = r*r*sintheta*sintheta;
}

void MonopoleGravity::MetricRadialDerivatives(
     Real r, Real theta, Real phi, Real Psi, Real dm, Real bh_mass,
     Real dPsi_dr, Real ddm_dr,
     Real &d1_g00, Real &d1_g11) const {
  
  const Real f = 1.0 - 2.0*bh_mass/r;
  Real r2 = SQR(r);
  Real d1_h_00 = -2.0*dm/r2 + 2.0/r*ddm_dr + 4.0*bh_mass/r2*Psi + 2.0*f*dPsi_dr;
  Real d1_h_11 = 2.0/(r*f*f)*(ddm_dr - (f+4.0*bh_mass/r)*dm/(r*f));
  
  Real d1_g_00 = -2.0*bh_mass / r2 + d1_h_00;
  Real d1_g_11 = -2.0*bh_mass / (r2*f*f) + d1_h_11;
  
  d1_g00 = d1_g_00;
  d1_g11 = d1_g_11;
}

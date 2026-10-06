
// C++ headers
#include <cmath>
#include <iostream>

//
#include "rgps_gravity.hpp"
#include "../athena.hpp"

void RGPSGravity::ConstructCovariantMetric(
    Real r, Real theta, Real phi, Real Phi, Real mgrav,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const {

  const Real sintheta = std::sin(theta);
  const Real alpha_sq = std::exp(2.0*Phi);
  const Real X_sq = 1/(1.0 - 2.0*mgrav/r);

  g00 = -alpha_sq;
  g01 = 0.0;
  g02 = 0.0;
  g03 = 0.0;
  
  g11 = X_sq;
  g12 = 0.0;
  g13 = 0.0;
  
  g22 = r*r;
  g23 = 0.0;
  g33 = r*r*sintheta*sintheta;
}

void RGPSGravity::MetricRadialDerivatives(
     Real r, Real theta, Real phi, Real Phi, Real mgrav,
     Real dPhi_dr, Real dmgrav_dr,
     Real &d1_g00, Real &d1_g11) const {
  
  const Real f = 1.0 - 2.0*mgrav/r;
  const Real r2= r*r;

  d1_g00 = std::exp(2.0*Phi) * 2.0*dPhi_dr;
  d1_g11 = -1.0/(f*f) * (2.0*mgrav/r2 - 2.0/r * dmgrav_dr);

}

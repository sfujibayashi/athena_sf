#ifndef INPUTS_OUTFLOW_BOUNDARY_DATA_HPP_
#define INPUTS_OUTFLOW_BOUNDARY_DATA_HPP_

#include <string>
#include <vector>

#include "../athena.hpp" // Real
#include "../utils/interp_table.hpp" // InterpTable2D

struct OutflowState {
  Real rho;
  Real press;
  Real vx;
  Real vy;
  Real vz;
  Real ye;
  Real entropy;
};

enum {
  IRHO_OUT = 0,
  IPRESS_OUT,
  IVX_OUT,
  IVY_OUT,
  IVZ_OUT,
  IYE_OUT,
  IENTROPY_OUT,
  NVAR_OUT
};

class OutflowBoundaryData {

public:
  OutflowBoundaryData(const std::string &filename);
  ~OutflowBoundaryData();
  
  OutflowState Interpolate(Real time, Real theta) const;

  Real GetTimeMin() const;
  Real GetTimeMax() const;
  Real GetThetaMin() const;
  Real GetThetaMax() const;

private:
  int ntime, ntheta;
  Real time_min, time_max;
  Real theta_min, theta_max;
  
  InterpTable2D table_;
};

#endif

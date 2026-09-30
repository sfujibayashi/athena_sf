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
  Real alpha;
  Real h;
  Real hut;
  Real psi;
  Real qe;
  Real rho_star;
  Real temp;
  Real ut;
  Real w;
};

enum {
  IRHO_OUT = 0,
  IPRESS_OUT,
  IVX_OUT,
  IVY_OUT,
  IVZ_OUT,
  IYE_OUT,
  IENTROPY_OUT,
  IALPHA_OUT,
  IENTHALPY_OUT,
  IHUT_OUT,
  IPSI_OUT,
  IQE_OUT,
  IRHOSTAR_OUT,
  IUT_OUT,
  ITEMP_OUT,
  IW_OUT,
  NVAR_OUT
};

class OutflowBoundaryData {

public:
  OutflowBoundaryData(const std::string &filename);
  ~OutflowBoundaryData();
  
  OutflowState Interpolate(Real time, Real theta) const;
  OutflowState GetState(int it, int j) const;
  
  void Analyze() const;

  Real GetTimeMin() const;
  Real GetTimeMax() const;
  Real GetThetaMin() const;
  Real GetThetaMax() const;

  Real GetNTheta() const;
  Real GetNTime() const;

private:
  int ntime, ntheta;
  Real time_min, time_max;
  Real theta_min, theta_max;

  Real r_ext_;
  Real h_min_global_;

  InterpTable2D table_;
};

#endif

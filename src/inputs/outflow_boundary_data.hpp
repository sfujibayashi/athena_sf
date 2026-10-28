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

struct AveragedOutflowState {
  // Used for primitives
  Real rho, press;       // cgs
  Real u1hat, u2hat, u3hat;  // local orthonormal frame u^(hat i)
  // Used for pscalar primitives
  Real ye;
  // Used for analysis
  Real E;
  Real alpha, psi;
};

enum{
  IRHO_AVE = 0,
  IPRESS_AVE,
  IU1HAT_AVE,
  IU2HAT_AVE,
  IU3HAT_AVE,
  IYE_AVE,
  IE_AVE,
  IALPHA_AVE,
  IPSI_AVE,
  NVAR_AVE
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

enum{
  IMASS_IN = 0,
  IMASS_OUT,
  IMEJ_GEOM,
  IMEJ_BERN,
  IMEJ_BIND,
  IEEJ_GEOM,
  IEEJ_BERN,
  IEEJ_BIND,
  IALPHA_AV,
  IALPHA_MIN,
  IALPHA_MAX,
  IPSI_AV,
  IPSI_MIN,
  IPSI_MAX,
  IMDOT_IN,
  IMDOT_OUT,
  IMDOT_GEOM,
  IMDOT_BERN,
  IMDOT_BIND,
  IEDOT_GEOM,
  IEDOT_BERN,
  IEDOT_BIND,
  IHNR_AV,
  IHATH_AV,
  NVAR_ANA
};


class OutflowBoundaryData {

public:
  OutflowBoundaryData(const std::string &filename);
  ~OutflowBoundaryData();
  
  OutflowState Interpolate(Real time, Real theta) const;
  OutflowState GetState(int it, int j) const;
  
  void Analyze(Real gamma_ad) const;

  Real GetTimeMin() const;
  Real GetTimeMax() const;
  Real GetThetaMin() const;
  Real GetThetaMax() const;

  int GetNTheta() const;
  int GetNTime() const;
  
  Real GetExtractionRadius() const;
  Real GetHMin() const;

  AveragedOutflowState GetAveragedData(int it, Real theta1, Real theta2) const;
  AveragedOutflowState GetPointData(int it, Real theta) const;
  Real GetTimeCGS(int it) const;

private:
  int ntime, ntheta;
  Real time_min, time_max;
  Real theta_min, theta_max;

  Real r_ext_;
  Real h_min_global_;

  InterpTableNonUniform2D table_;
  AthenaArray<Real> time_, theta_c_, theta_d_, theta_u_;
};



class OutflowBoundaryTable {
public:
  OutflowBoundaryTable(const OutflowBoundaryData &source,
                       const AthenaArray<Real> &theta_faces,
		       const int js, const int je);
  
  void Interpolate(Real time_cgs, int js, int je, AthenaArray<Real> &out) const;
  AveragedOutflowState InterpolatePoint(Real time_cgs, int j) const;
  Real GetTime(const int it) const;
  Real GetTimeMax() const;
  int GetNTime() const;
  void Analyze(Real gamma_ad,
	       const AthenaArray<Real> &theta,
	       const AthenaArray<Real> &theta_faces,
	       const int js, const int je,
	       AthenaArray<Real> &out) const;
  
private:
  int ntime_;
  AthenaArray<Real> time_;
  AthenaArray<Real> data_;

  Real r_ext_;
  Real h_min_global_;
};

#endif

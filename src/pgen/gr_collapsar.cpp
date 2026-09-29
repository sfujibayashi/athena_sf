//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast_helmholtz.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// progenitor-reader
#include "../inputs/progenitor_reader.hpp"
#include "../metric/metric.hpp"
#include "../gravity/monopole_gravity.hpp"
#include "../inputs/outflow_boundary_data.hpp"

struct CollapsedProfile {
  int nface = 0;
  int ncell = 0;

  std::vector<Real> mass_face;
  std::vector<Real> radius_face;
  std::vector<Real> ur_face;

  std::vector<Real> mass;
  std::vector<Real> radius;

  std::vector<Real> rho;
  std::vector<Real> press;
  std::vector<Real> ur;

  std::vector<Real> jrot;

  Real EnclosedMass(Real r){
    int ilo, ihi;
    ilo = 0;
    ihi = nface-1;

    if (r>=radius_face[ihi]){
      return mass_face[ihi];
    }
    
    while(ihi-ilo>=2){
      int i=(ihi+ilo)/2;
      if ( radius_face[i] <= r ){
        ilo = i;
      }else{
        ihi = i;
      }
    }
    return mass_face[ilo]; // not interpolated for now.
  }

  int FirstCellIndex(){
    int i = 0;
    while(radius[i]==0.0){
      i++;
    }
    return i;
  }

  int CellIndexFromRadius(Real r){
    int ilo, ihi;
    ihi = ncell-1;

    ilo = FirstCellIndex();
    
    if (r>=radius_face[nface-1]){
      return ncell;
    } else if ( r>= radius[ncell-1]){
      return ncell-1;
    } else if ( r<radius[ilo]){
      return -1;
    }
    
    while(ihi-ilo>=2){
      int i=(ihi+ilo)/2;
      if ( radius[i] <= r ){
        ilo = i;
      }else{
        ihi = i;
      }
    }

    return ilo;
    
  }
};

Real GetTimeFromBlackHoleMass(const ProgenitorProfile &progenitor, Real bh_mass);
CollapsedProfile CollapseProgenitor(const ProgenitorProfile &progenitor, Real t0);

void InjectionInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                      Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  CollapsedProfile collapsed;  
  Real gamma_gas;

  Real timescale_cut;
  Real W_max;

  // Real HistoryBlackHoleMass(MeshBlock *pmb, int iout);
  // Real HistoryBlackHoleMassAccretionRate(MeshBlock *pmb, int iout);
  
  Real HistoryBlackHoleMass(MeshBlock *pmb, int iout) {
    const Real mass_to_length =
      pmb->pmy_mesh->punit->grav_const_code
      / SQR(pmb->pmy_mesh->punit->speed_of_light_code);
    
    return pmb->pmetric->GetBlackHoleMass()/mass_to_length;
  }
  
  Real HistoryBlackHoleMassAccretionRate(MeshBlock *pmb, int iout) {
    return pmb->pmy_mesh->pmonograv->mdot_bh_;
  }

  Real HistoryOuterMassFlux(MeshBlock *pmb, int iout) {
    // This block does not touch the physical outer-x1 boundary.
    if (pmb->pbval->block_bcs[BoundaryFace::outer_x1]
        == BoundaryFlag::block) {
      return 0.0;
    }
    
    Real mdot_out = 0.0;
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        const Real area =
          pmb->pcoord->GetFace1Area(k, j, pmb->ie+1);
        
        mdot_out += area
          * pmb->phydro->flux[X1DIR](IDN,k,j,pmb->ie+1);
      }
    }
    
    return mdot_out;
  }


  Real HistorySumAbs2Mom(MeshBlock *pmb, int iout) {
    // This block does not touch the physical outer-x1 boundary.
    Real sum_2mom = 0.0;
    
    AthenaArray<Real> vol;
    vol.NewAthenaArray(pmb->ie-pmb->is+1);

    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          sum_2mom += vol(i) * std::abs(pmb->phydro->u(IM2,k,j,i));
        }
      }
    }
    
    return sum_2mom;
  }

}

void Mesh::InitUserMeshData(ParameterInput *pin) {

  std::string filename = pin->GetString("problem", "progenitor_file");

  const bool initialize_with_mass =
    pin->GetOrAddBoolean("coord", "initialize_with_mass", true);

  ProgenitorProfile progenitor = ReadProgenitorProfile(filename);

  Real t0=0.0;
  if(initialize_with_mass){
    Real bh_mass_init = pin->GetReal("problem", "initial_bh_mass");
    t0 = GetTimeFromBlackHoleMass(progenitor, bh_mass_init);
  }else{
    t0 = pin->GetReal("problem", "t0");
  }
  gamma_gas = pin->GetReal("hydro", "gamma");
  
  
  collapsed = CollapseProgenitor(progenitor, t0);
  
  Real rin_code = mesh_size.x1min;
  Real rin_cgs = rin_code*punit->code_length_cgs;
  Real M_inner_cgs = collapsed.EnclosedMass(rin_cgs);
  Real J_inner_cgs = 0.0;

  Real m_bh_code = Constants::grav_const_cgs * M_inner_cgs
    / SQR(Constants::speed_of_light_cgs)
    / punit->code_length_cgs;

  Real ang_bh_code = J_inner_cgs
    * std::pow(Constants::grav_const_cgs / SQR(Constants::speed_of_light_cgs)
    / punit->code_length_cgs, 2);
  
  printf("Black hole mass (cgs,code unit)=%15.7e, %15.7e\n",M_inner_cgs, m_bh_code);
  // pin->SetReal("coord", "m", m_bh_code);
  // pin->SetReal("coord", "j", ang_bh_code);
  // m_bh_code = pin->GetReal("coord", "m");

  // std::cout << std::setprecision(17)
  //           << "m_bh_code = " << m_bh_code
  //           << "  m_from_pin = " << pin->GetReal("coord","m")
  //           << std::endl;

  // inject BC
  std::string fname =
    pin->GetString("problem", "outflow_file");
  poutflow = new OutflowBoundaryData(fname);

  timescale_cut = pin->GetReal("problem", "timescale_cut");
  W_max = pin->GetReal("problem", "W_max");
  
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1,
                             InjectionInnerX1);
  
  // output
  AllocateUserHistoryOutput(3);
  
  EnrollUserHistoryOutput(0, HistoryBlackHoleMass, "m_bh",
                          UserHistoryOperation::max);
  
  EnrollUserHistoryOutput(1, HistoryBlackHoleMassAccretionRate, "mdot_bh",
                          UserHistoryOperation::max);

  EnrollUserHistoryOutput(2, HistoryOuterMassFlux, "mdot_out",
                          UserHistoryOperation::sum);

  // EnrollUserHistoryOutput(3, HistorySumAbs2Mom, "Abs 2-Mom",
  //                         UserHistoryOperation::sum);

  //
  AllocateRealUserMeshDataField(2);
  // BH mass
  ruser_mesh_data[0].NewAthenaArray(1);
  // BH spin
  ruser_mesh_data[1].NewAthenaArray(1);

  // canonical BH mass
  ruser_mesh_data[0](0) = m_bh_code;
  // canonical BH spin
  ruser_mesh_data[1](0) = ang_bh_code;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(9);
  
  SetUserOutputVariableName(0, "gtt");
  SetUserOutputVariableName(1, "grr");
  SetUserOutputVariableName(2, "alpha-1");
  SetUserOutputVariableName(3, "Lorentz-1");
  SetUserOutputVariableName(4, "u_t+1");
  SetUserOutputVariableName(5, "enthalpy-1");
  SetUserOutputVariableName(6, "delta_m");
  SetUserOutputVariableName(7, "Phi");
  SetUserOutputVariableName(8, "q");

  //
  AllocateRealUserMeshBlockDataField(2);
  
  // for Psi_face1
  ruser_meshblock_data[0].NewAthenaArray(ncells1 + 1);
  // for delta_m_face1
  ruser_meshblock_data[1].NewAthenaArray(ncells1 + 1);

  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();
}

void MeshBlock::UserWorkInLoop(void) {

}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie + NGHOST + 1);
  gi.NewAthenaArray(NMETRIC, ie + NGHOST + 1);

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);

      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = g(I00,i);
        user_out_var(1,k,j,i) = g(I11,i);
        user_out_var(2,k,j,i) = std::sqrt(-1.0 / gi(I00,i)) - 1.0;
      }
    }
  }

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->CellMetric(k, j, is, ie, g, gi);

      for (int i = is; i <= ie; ++i) {
        Real uu1 = phydro->w(IVX,k,j,i);
        Real uu2 = phydro->w(IVY,k,j,i);
        Real uu3 = phydro->w(IVZ,k,j,i);
        
        Real tmp = g(I11,i)*SQR(uu1)
          + 2.0*g(I12,i)*uu1*uu2
          + 2.0*g(I13,i)*uu1*uu3
          + g(I22,i)*SQR(uu2)
          + 2.0*g(I23,i)*uu2*uu3
          + g(I33,i)*SQR(uu3);
        Real lorentz = std::sqrt(1.0 + tmp);
        user_out_var(3,k,j,i) = lorentz-1.0;

        Real alpha = std::sqrt(-1.0 / gi(I00,i));
        Real u_t = -alpha*lorentz;
        user_out_var(4,k,j,i) = u_t+1.0;

        Real press = phydro->w(IPR,k,j,i);
        Real rho   = phydro->w(IDN,k,j,i);
        Real enthalpy = 1.0 + gamma_gas/(gamma_gas - 1.0) * press/rho;
        user_out_var(5,k,j,i) = enthalpy-1.0;
      }
    }
  }


  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        user_out_var(6,k,j,i) = pmetric->CellDeltaM(i);
        user_out_var(7,k,j,i) = pmetric->CellPsi(i);
        user_out_var(8,k,j,i) = pmetric->CellDensitizationFactor(k,j,i);
        
      }
    }
  }
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") != 0 && 
      std::strcmp(COORDINATE_SYSTEM, "minkowski") != 0 && 
      std::strcmp(COORDINATE_SYSTEM, "gr_dynamic") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_collapsar.cpp" << std::endl
        << "This problem generator requires Schwarzschild/Minkowski coordinates."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  
  Real rho_atmos = pin->GetReal("problem", "rho_atmos");
  Real press_atmos = pin->GetReal("problem", "press_atmos");
  int ind_first = collapsed.FirstCellIndex();

  // Prepare scratch arrays
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie+NGHOST+1);
  gi.NewAthenaArray(NMETRIC, ie+NGHOST+1);
  
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i=is; i<=ie; i++) {
        Real rad_cgs = pcoord->x1v(i)*pmy_mesh->punit->code_length_cgs;
        
        int ind = collapsed.CellIndexFromRadius(rad_cgs);
        
        Real rho_cgs, press_cgs, uu1, uu2, uu3;
        uu2 = 0.0;
        uu3 = 0.0;
        if( ind>=collapsed.ncell ){
          // outside the star atmosphere
          rho_cgs = rho_atmos;
          press_cgs = press_atmos;
          uu1 = 0.0;
        }else if (ind==collapsed.ncell-1){
          // between the last cell center and stellar surface
          rho_cgs   = collapsed.rho[ind];
          press_cgs = collapsed.press[ind];
          uu1 = collapsed.ur[ind] / Constants::speed_of_light_cgs;
        }else if (ind<0){
          rho_cgs = collapsed.rho[ind_first];
          press_cgs= collapsed.press[ind_first];
          uu1 = collapsed.ur[ind_first]/Constants::speed_of_light_cgs;
        } else {
          Real xx1 = (rad_cgs-collapsed.radius[ind])/(collapsed.radius[ind+1]-collapsed.radius[ind]);
          Real xx0 = 1.0-xx1;
          rho_cgs  = xx0*collapsed.rho[ind] + xx1*collapsed.rho[ind+1];
          press_cgs= xx0*collapsed.press[ind] + xx1*collapsed.press[ind+1];
          uu1 = (xx0*collapsed.ur[ind] + xx1*collapsed.ur[ind+1])/Constants::speed_of_light_cgs;
        }
        Real rho_code = rho_cgs/pmy_mesh->punit->code_density_cgs;
        Real press_code = press_cgs/pmy_mesh->punit->code_pressure_cgs;
        
        phydro->w(IDN,k,j,i) = rho_code;
        phydro->w(IPR,k,j,i) = press_code;
        phydro->w(IVX,k,j,i) = uu1;
        phydro->w(IVY,k,j,i) = uu2;
        phydro->w(IVZ,k,j,i) = uu3;

        // Initial guess / previous primitives
        phydro->w1(IDN,k,j,i) = rho_code;
        phydro->w1(IPR,k,j,i) = press_code;
        phydro->w1(IVX,k,j,i) = uu1;
        phydro->w1(IVY,k,j,i) = uu2;
        phydro->w1(IVZ,k,j,i) = uu3;
        printf("i, rho, press, uu1 = %5d %12.4e %12.4e %12.4e\n", i, rho_code, press_code, uu1);
      }
    }
  }

  const bool initial_metric_update =
    pin->GetOrAddBoolean("coord", "initial_metric_update", true);
  
  const bool time_metric_update =
    pin->GetOrAddBoolean("coord", "time_metric_update", true);
  
  if (time_metric_update && !initial_metric_update) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_collapsar.cpp" << std::endl
        << "time_metric_update=true requires initial_metric_update=true."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  
  // Convert primitive -> conserved
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  bb.ZeroClear();

  peos->PrimitiveToConserved(
      phydro->w, bb, phydro->u, pcoord,
      is, ie, js, je, ks, ke);

}


// return eta that satisfies
// t_collapse = t - t_m0 = sqrt(r0^3/8Gm)*(eta + sin(eta)).
// r(t) = (1/2)*r0*(1+cos(eta))
// eta = 0 if t - t_m0 < 0.
Real GetCollapsedRadius(const Real t_collapse, const Real r0, const Real m){
  Real tscale = std::sqrt(std::pow(r0,3)/(8.0*Constants::grav_const_cgs*m));
  Real tff= PI*tscale;
  Real eta_lo = 0.0;
  Real eta_hi = PI;

  Real tol = 1.0e-10;

  if(t_collapse < 0.0){
    return r0;
  }else if(t_collapse > tff){
    return 0.0;
  }

  Real t_tilde = t_collapse/tscale;
  Real f_hi = eta_hi + std::sin(eta_hi) - t_tilde;
  Real f_lo = eta_lo + std::sin(eta_lo) - t_tilde;

  while( std::abs(eta_lo - eta_hi)>tol){
    Real eta = 0.5*(eta_lo+eta_hi);
    Real f = eta + std::sin(eta) - t_tilde;
    if( f_hi*f >= 0.0){
      f_hi = f;
      eta_hi = eta;
    }else{
      f_lo = f;
      eta_lo = eta;
    }
  }
  
  Real eta= 0.5*(eta_hi + eta_lo);
  return r0*0.5*(1.0+std::cos(eta));
}

Real GetTimeFromBlackHoleMass(const ProgenitorProfile &progenitor, Real bh_mass){

  const Real Msun = Constants::solar_mass_cgs;
  const Real G = Constants::grav_const_cgs;
  const Real c = Constants::speed_of_light_cgs;
  
  // find radius for m=bh_mass  
  const int nface = progenitor.nface;
  int ilo, ihi;
  ilo = 0;
  ihi = nface-1;
  
  if (bh_mass>=progenitor.mass_face[ihi]/Msun){
    std::cout << "bh_mass is larger than the total mass of the star."<< std::endl;
    return -1.0;
  }
  
  while(ihi-ilo>=2){
    int i=(ihi+ilo)/2;
    if ( progenitor.mass_face[i]/Msun <= bh_mass ){
      ilo = i;
    }else{
      ihi = i;
    }
  }
  // Radius of enclosed mass = bh_mass
  Real r_m0 = progenitor.radius_face[ilo]; // not interpolated for now.

  // derive collapse parameter eta for which r(m) = 2*m.
  Real eta = std::acos(4.0*G*bh_mass*Msun/(c*c) - 1.0);

  // derive sound crossing time
  Real t_m0 = 0.0;
  for(int i=0; i<ilo; ++i){
    Real dr = progenitor.radius_face[i+1] - progenitor.radius_face[i];
    Real dt = dr/progenitor.csound[i];
    t_m0 += dt;
  }

  Real tau = t_m0 + std::sqrt(r_m0*r_m0*r_m0/(8.0*G*bh_mass*Msun)) * (eta + std::sin(eta));

  return tau;
}

CollapsedProfile CollapseProgenitor(const ProgenitorProfile &progenitor, Real t0){
  CollapsedProfile collapsed;
  std::vector<Real> t_m0;

  collapsed.ncell = progenitor.ncell;
  collapsed.nface = progenitor.nface;

  t_m0.resize(collapsed.nface);

  collapsed.mass_face = progenitor.mass_face;
  
  collapsed.radius_face.resize(collapsed.nface);
  collapsed.ur_face.resize(collapsed.nface);

  collapsed.radius.resize(collapsed.ncell);
  collapsed.rho.resize(collapsed.ncell);
  collapsed.press.resize(collapsed.ncell);
  collapsed.ur.resize(collapsed.ncell);
  
  collapsed.mass = progenitor.mass;
  collapsed.jrot = progenitor.jrot;


  t_m0[0] = 0.0;

  for(int i=0; i<collapsed.ncell; ++i){
    Real dr = progenitor.radius_face[i+1] - progenitor.radius_face[i];
    Real dt = dr/progenitor.csound[i];
    t_m0[i+1] = t_m0[i] + dt;
  }

  collapsed.radius_face[0] = 0.0;
  collapsed.ur_face[0] = 0.0;
  
  for (int i=1; i<collapsed.nface; ++i) {
    Real t_collapse = t0 - t_m0[i];
    
    if (t_collapse <= 0.0) {
      // sound wave has not reached this shell
      collapsed.radius_face[i] = progenitor.radius_face[i];
      collapsed.ur_face[i] = 0.0;
    } else {
      Real rf0 = progenitor.radius_face[i];
      Real mf  = collapsed.mass_face[i];
      Real rf = GetCollapsedRadius(t_collapse, rf0, mf);
      collapsed.radius_face[i] = rf;
      if(rf>0.0){
        collapsed.ur_face[i] = -std::sqrt(2.0*Constants::grav_const_cgs*mf*(rf0 - rf)/(rf0*rf));
      }else{
        collapsed.ur_face[i] = 0.0;
      }
      
    }
  }
  
  for (int i=0; i<collapsed.ncell; ++i){
    Real rp = collapsed.radius_face[i+1];
    Real rm = collapsed.radius_face[i];
    
    collapsed.radius[i] = std::cbrt(0.5*(rm*rm*rm + rp*rp*rp));
    collapsed.ur[i] = 0.5*(collapsed.ur_face[i] + collapsed.ur_face[i+1]);
    
    Real dv = 4.0*PI/3.0*(rp*rp*rp - rm*rm*rm);
    Real dm = collapsed.mass_face[i+1] - collapsed.mass_face[i];

    if (dv > 0.0){
      collapsed.rho[i] = dm/dv;
      collapsed.press[i] = progenitor.press[i]*std::pow(collapsed.rho[i]/progenitor.rho[i], gamma_gas);
    }else{
      collapsed.ur[i] = 0.0;
      collapsed.rho[i] = 0.0;
      collapsed.press[i] = 0.0;
    }
  }
  
  return collapsed;
}

void InjectionInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                      Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ku, int ngh){

  const Real code_time_cgs = pmb->pmy_mesh->punit->code_time_cgs;
  const Real code_rho_cgs = pmb->pmy_mesh->punit->code_density_cgs;
  const Real code_press_cgs = pmb->pmy_mesh->punit->code_pressure_cgs;

  const Real vv_max = 1.0 - 1.0/(W_max*W_max);
  Real time_cgs = time * code_time_cgs;

  Real fac = 1.0;
  Real time_interp = time_cgs;

  // Avoid extrapolation
  // should determine how we will deal with after running out of data.
  if      ( time_cgs > pmb->pmy_mesh->poutflow->GetTimeMax() ){
    time_interp = pmb->pmy_mesh->poutflow->GetTimeMax();
    fac = std::exp(- (time_cgs-pmb->pmy_mesh->poutflow->GetTimeMax())/timescale_cut );
  }else if( time_cgs < pmb->pmy_mesh->poutflow->GetTimeMin() ){
    time_interp = pmb->pmy_mesh->poutflow->GetTimeMin();
    fac = 1.0;
  }
  
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    // Real phi = pmb->pcoord->x3v(k);
    // Real cosphi = std::cos(phi);
    // Real sinphi = std::sin(phi);
    for (int j=jl; j<=ju; ++j) {
      Real theta = pmb->pcoord->x2v(j);

      // Avoid extrapolation.
      // We rotate the data with theta_max/min assuming v^r,theta,phi do not depend on theta
      // in angle range theta_max < theta < pi/2.
      Real theta_interp = theta;
      if      (theta < pmb->pmy_mesh->poutflow->GetThetaMin()){
        theta_interp = pmb->pmy_mesh->poutflow->GetThetaMin();
      }else if(theta > pmb->pmy_mesh->poutflow->GetThetaMax()){
        theta_interp = pmb->pmy_mesh->poutflow->GetThetaMax();
      }

      Real costheta = std::cos(theta_interp);
      Real sintheta = std::sin(theta_interp);
      // Real sin2theta = sintheta*sintheta;
      
      OutflowState state = pmb->pmy_mesh->poutflow->Interpolate(time_interp, theta_interp);

      // // Cartesian to spherical-polar matrix.
      // Real m1_x = sintheta*cosphi;
      // Real m1_y = sintheta*sinphi;
      // Real m1_z = costheta;
      
      // // should be divided by r additionally.
      // Real m2_x = costheta*cosphi;
      // Real m2_y = costheta*sinphi;
      // Real m2_z = -sintheta;

      // // should be divided by r additionally.
      // Real m3_x = -sinphi/sintheta;
      // Real m3_y = cosphi/sintheta;
      // Real m3_z = 0.0;

      
      // Cylindrical (y=0) to spherical-polar matrix.
      Real m1_R = sintheta;
      Real m1_z = costheta;
      
      // should be divided by r additionally.
      Real m2_R = costheta;
      Real m2_z = -sintheta;

      // should be divided by r additionally.
      Real m3_y = 1.0/sintheta;

#pragma omp simd
      for (int n=1; n<=ngh; ++n) {
        int i = il - n;
        Real r = pmb->pcoord->x1v(i);
        // Real r2=r*r;
        
        // Construct metric
        Real g00, g01, g02, g03;
        Real g11, g12, g13, g22, g23, g33;
        pmb->pmetric->ConstructCellCovariantMetric(k,j,i,
                                          g00, g01, g02, g03,
                                          g11, g12, g13, g22, g23, g33);
        Real alpha = std::sqrt(-g00);

        // assume beta^i = 0.
        Real v1 = (m1_R*state.vx + m1_z*state.vz)   / alpha;
        Real v2 = (m2_R*state.vx + m2_z*state.vz)/r / alpha;
        Real v3 = (m3_y*state.vy)/r / alpha;
        
        Real vv = g11*v1*v1 + g22*v2*v2 + g33*v3*v3
           + 2.0*(g12*v1*v2 + g13*v1*v3 + g23*v2*v3);

        if (vv > vv_max){
          std::cout << "Warning: v^2 > v2_max." << vv
                    << " at t=" << time_interp << "theta=" << theta_interp << std::endl;
          const Real factor = std::sqrt(vv_max/vv);
          v1 *= factor;
          v2 *= factor;
          v3 *= factor;
          vv = vv_max;
        }
        // lorentz factor
        Real W = 1.0/std::sqrt(1.0-vv);

        // Athena's primitive.
        Real uu1 = W * v1;
        Real uu2 = W * v2;
        Real uu3 = W * v3;
        
        prim(IDN,k,j,i) = state.rho/code_rho_cgs * fac;
        prim(IPR,k,j,i) = state.press/code_press_cgs * fac;
        prim(IVX,k,j,i) = uu1;
        prim(IVY,k,j,i) = uu2;
        prim(IVZ,k,j,i) = uu3;
      }
    }
  }
}

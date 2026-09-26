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

CollapsedProfile CollapseProgenitor(const ProgenitorProfile &progenitor, Real t0);

namespace {
  Real m_bh, a_bh;
  CollapsedProfile collapsed;  
  Real gamma_gas;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {

  std::string filename = pin->GetString("problem", "progenitor_file");
  Real t0 = pin->GetReal("problem", "t0");

  gamma_gas = pin->GetReal("hydro", "gamma");
  
  ProgenitorProfile progenitor = ReadProgenitorProfile(filename);
  
  collapsed = CollapseProgenitor(progenitor, t0);
  
  Real rin_code = mesh_size.x1min;
  Real rin_cgs = rin_code*punit->code_length_cgs;
  Real M_inner_cgs = collapsed.EnclosedMass(rin_cgs);

  Real m_bh_code = Constants::grav_const_cgs * M_inner_cgs
    / SQR(Constants::speed_of_light_cgs)
    / punit->code_length_cgs;
  printf("Black hole mass (cgs,code unit)=%15.7e, %15.7e\n",M_inner_cgs, m_bh_code);
  pin->SetReal("coord", "m", m_bh_code);
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
        user_out_var(6,k,j,i) = pmetric->delta_m_(i);
        user_out_var(7,k,j,i) = pmetric->Psi_(i);
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
  // Calculate metric perturbation from primitive rho.
  pmetric->Update();

  std::cout << "r=" << pcoord->x1v(is)
            << " dm=" << pmetric->delta_m_(is)
            << " Psi=" << pmetric->Psi_(is)
            << " q=" << pmetric->CellDensitizationFactor(0,0,is)
            << std::endl;

  // Convert primitive -> conserved
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  bb.ZeroClear();

  peos->PrimitiveToConserved(
      phydro->w, bb, phydro->u, pcoord,
      is, ie, js, je, ks, ke);
  //std::abort();
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

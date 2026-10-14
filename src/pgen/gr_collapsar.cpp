//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_collapsar.cpp
//! \brief Problem generator for long-term simulation of collapse and explosion of stars.

// python3 configure.py --prob gr_collapsar --coord gr_dynamic --eos adiabatic --flux hllc --nghost 2 --nscalars 22 -g -t -mpi -hdf5 --hdf5_path /home/sfujibayashi/libs/hdf5/2.0.0 --cxx g++ 

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
#include "../gravity/dynamic_metric_driver.hpp"
#include "../inputs/outflow_boundary_data.hpp"

#include "gr_collapsar_history.hpp"

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

  
  std::vector<Real> ye;
  std::array<std::vector<Real>, ProgenitorSpecies::NPROG_SPECIES> x;
  
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

    Real mass;
    if (radius_face[ilo] <= 0.0){
      Real w1 = (r-radius_face[ilo])/(radius_face[ihi] - radius_face[ilo]);
      Real w0 = 1.0 - w1;
      mass = w0*mass_face[ilo] + w1*mass_face[ihi]; // linear-linear interpolation
    }else{
      Real w1 = (std::log(r)-std::log(radius_face[ilo]))/(std::log(radius_face[ihi]) - std::log(radius_face[ilo]));
      Real w0 = 1.0 - w1;
      mass = std::exp(w0*std::log(mass_face[ilo]) + w1*std::log(mass_face[ihi])); //log-log interpolation
    }

    return mass;
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
// void DiodeInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                       Real time, Real dt,
//                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                       Real time, Real dt,
//                       int il, int iu, int jl, int ju, int kl, int ku, int ngh);

Real LogarithmicX1(Real x, RegionSize rs);

namespace {
  CollapsedProfile collapsed;  
  Real gamma_gas;

  Real timescale_cut;
  Real W_max;

  
  enum ScalarIndex {
    IYE = 0,
    IXINJ,
    IXPHYS,
    
    IXNEUT,
    IXH1,
    IXHE3,
    IXHE4,
    IXC12,
    IXN14,
    IXO16,
    IXNE20,
    IXMG24,
    IXSI28,
    IXS32,
    IXAR36,
    IXCA40,
    IXTI44,
    IXCR48,
    IXCR56,
    IXFE52,
    IXFE54,
    IXFE56,
    IXNI56,

    
    NSCALAR_REQUIRED
  };
  
  const int prog_to_scalar[ProgenitorSpecies::NPROG_SPECIES] = {
    IXNEUT,  // IPROG_NEUT
    IXH1,    // IPROG_PROT
    IXH1,    // IPROG_H1
    IXHE3,   // IPROG_HE3
    IXHE4,   // IPROG_HE4
    IXC12,   // IPROG_C12
    IXN14,   // IPROG_N14
    IXO16,   // IPROG_O16
    IXNE20,  // IPROG_NE20
    IXMG24,  // IPROG_MG24
    IXSI28,  // IPROG_SI28
    IXS32,   // IPROG_S32
    IXAR36,  // IPROG_AR36
    IXCA40,  // IPROG_CA40
    IXTI44,  // IPROG_TI44
    IXCR48,  // IPROG_CR48
    IXCR56,  // IPROG_CR56
    IXFE52,  // IPROG_FE52
    IXFE54,  // IPROG_FE54
    IXFE56,  // IPROG_FE56
    IXNI56   // IPROG_NI56
  };
  
  // Real HistoryBlackHoleMass(MeshBlock *pmb, int iout);
  // Real HistoryBlackHoleMassAccretionRate(MeshBlock *pmb, int iout);

  Real CellSpecificBindingEnergy(MeshBlock *pmb, int k, int j, int i,
			     const AthenaArray<Real> &gi) {

#if VALENCIA_FORMULATION

    const Real qD   = pmb->phydro->u(IDN,k,j,i);
    const Real qtau = pmb->phydro->u(IEN,k,j,i);
    const Real qS_1  = pmb->phydro->u(IM1,k,j,i);
    const Real qS_2  = pmb->phydro->u(IM2,k,j,i);
    const Real qS_3  = pmb->phydro->u(IM3,k,j,i);

    const Real alpha = std::sqrt(-1.0/gi(I00,i));

    const Real beta1 = alpha*alpha*gi(I01,i);
    const Real beta2 = alpha*alpha*gi(I02,i);
    const Real beta3 = alpha*alpha*gi(I03,i);

    return (alpha*(qD + qtau)
	    - beta1*qS_1 - beta2*qS_2 - beta3*qS_3)/qD - 1.0;

#else

    const Real rho_ut = pmb->phydro->u(IDN,k,j,i);
    const Real Tt_t   = pmb->phydro->u(IEN,k,j,i);

    return -Tt_t/rho_ut - 1.0;

#endif
  }

  
  Real HistoryBlackHoleMass(MeshBlock *pmb, int iout) {
    const Real mass_to_length =
      pmb->pmy_mesh->punit->grav_const_code
      / SQR(pmb->pmy_mesh->punit->speed_of_light_code);
    
    return pmb->pmy_mesh->pmetric_driver->GetBlackHoleMass()
      / mass_to_length / pmb->pmy_mesh->punit->solar_mass_code;
  }
  
  Real HistoryBlackHoleMassAccretionRate(MeshBlock *pmb, int iout) {
    const Real mdot_geom =
      pmb->pmy_mesh->pmetric_driver->GetBlackHoleMassAccretionRate();
    const Real mass_to_length =
      pmb->pmy_mesh->punit->grav_const_code
      / SQR(pmb->pmy_mesh->punit->speed_of_light_code);

    return mdot_geom / mass_to_length
      / pmb->pmy_mesh->punit->solar_mass_code
      / pmb->pmy_mesh->punit->code_time_cgs;

  }

  Real HistoryMInnerBlocks(MeshBlock *pmb, int iout) {
    // This block does not touch the physical outer-x1 boundary.
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::block) {
      return 0.0;
    }

    AthenaArray<Real> vol;
    vol.NewAthenaArray(pmb->ncells1);
    Real mass = 0.0;
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
	for (int i=pmb->is; i<=pmb->ie; ++i) {
	  const Real qD = pmb->phydro->u(IDN,k,j,i);
	  mass += vol(i) * qD;
	}
      }
    }
    
    return mass / pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEInnerBlocks(MeshBlock *pmb, int iout) {
    // This block does not touch the physical outer-x1 boundary.
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::block) {
      return 0.0;
    }

    AthenaArray<Real> vol;
    vol.NewAthenaArray(pmb->ncells1);

    Real E = 0.0;
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);

	for (int i=pmb->is; i<=pmb->ie; ++i) {
#if VALENCIA_FORMULATION
	  const Real qD   = pmb->phydro->u(IDN,k,j,i);
	  const Real qtau = pmb->phydro->u(IEN,k,j,i);
	  const Real q = pmb->pmy_mesh->pmetric_driver->CellSpatialDensitizationFactor(pmb, k,j,i);
	  const Real r = pmb->pcoord->x1v(i);
	  const Real Eatmos = pmb->peos->GetEnergyFloor(r);
	  E += vol(i) * ((qD+qtau)/q-Eatmos);
#else
	  const Real Tt_t = pmb->phydro->u(IEN,k,j,i);
	  E += vol(i) * (-Tt_t);
#endif
	}
      }
    }
    
    return E * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      / pmb->pmy_mesh->punit->erg_code;
  }


  Real HistoryJInnerBlocks(MeshBlock *pmb, int iout) {
    // This block does not touch the physical outer-x1 boundary.
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::block) {
      return 0.0;
    }

    AthenaArray<Real> vol;
    vol.NewAthenaArray(pmb->ncells1);
    
    Real J = 0.0;
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
	for (int i=pmb->is; i<=pmb->ie; ++i) {
#if VALENCIA_FORMULATION
	  const Real qS_3 = pmb->phydro->u(IM3,k,j,i);
	  J += vol(i) * qS_3;
#else
	  const Real Tt_3 = pmb->phydro->u(IM3,k,j,i);
	  J += vol(i) * Tt_3;
#endif
	}
      }
    }
    
    return J * pmb->pmy_mesh->punit->code_velocity_cgs * pmb->pmy_mesh->punit->code_length_cgs;
  }


  enum EjectaHistoryCache {
    HCYCLE = 0,
    
    MGEOM,
    EGEOM,
    MBERN,
    EBERN,
    MBIND,
    EBIND,

    RMAX,
    
    MGEOM_OUT,
    MBAL,
    
    MTEXP,
    MTEXP2,
    
    MR2,
    MRVR,
    MVR2,

    NEJECTA_CACHE
  };

  enum BoundaryFluxCache {
    MDOT_IN_BC = 0,
    EDOT_IN_BC,
    MDOT_IN_NUM,
    EDOT_IN_NUM,
    MDOT_OUT_BC,
    EDOT_OUT_BC,
    MDOT_OUT_NUM,
    EDOT_OUT_NUM,

    BDOT_IN_BC,
    BDOT_IN_NUM,
    BDOT_OUT_BC,
    BDOT_OUT_NUM,

    MDOT_OUT,

    MDOT_INJ_IN,
    MDOT_INJ_OUT,
    
    NBOUNDARY_FLUX_CACHE
  };

  // enum BoundaryFluxIntegral {
  //   MASS_IN_BC = 0,
  //   ENERGY_IN_BC,
  //   MASS_IN_NUM,
  //   ENERGY_IN_NUM,
  //   MASS_OUT_BC,
  //   ENERGY_OUT_BC,
  //   MASS_OUT_NUM,
  //   ENERGY_OUT_NUM,

  //   MASS_OUT,

  //   MASS_REMOVED,
    
  //   NBOUNDARY_FLUX_INTEGRAL
  // };
  using namespace GRCollapsarHistory;
  
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

  void UpdateEjectaHistoryCache(MeshBlock *pmb) {
    AthenaArray<Real> &cache = pmb->ruser_meshblock_data[2];
    
    const int ncycle = pmb->pmy_mesh->ncycle;
    
    // Already evaluated at this timestep
    if (cache(HCYCLE) == static_cast<Real>(ncycle)) {
      return;
    }
    
    for (int n=1; n<NEJECTA_CACHE; ++n) {
      cache(n) = 0.0;
    }
    
    AthenaArray<Real> vol;
    vol.NewAthenaArray(pmb->ncells1);
    
    AthenaArray<Real> g, gi;
    g.NewAthenaArray(NMETRIC, pmb->ncells1);
    gi.NewAthenaArray(NMETRIC, pmb->ncells1);
    
    // Example; make this an input parameter later
    const Real chiP_cut = 1.0e-2;
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
	
	pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
	
	pmb->pcoord->CellMetric(k, j, pmb->is, pmb->ie, g, gi);

	for (int i=pmb->is; i<=pmb->ie; ++i) {

	  const Real rho = pmb->phydro->w(IDN,k,j,i);
	  const Real press = pmb->phydro->w(IPR,k,j,i);

	  const Real uu1 = pmb->phydro->w(IVX,k,j,i);
	  const Real uu2 = pmb->phydro->w(IVY,k,j,i);
	  const Real uu3 = pmb->phydro->w(IVZ,k,j,i);

	  const Real usq =
	    g(I11,i)*uu1*uu1
	    + 2.0*g(I12,i)*uu1*uu2
	    + 2.0*g(I13,i)*uu1*uu3
	    + g(I22,i)*uu2*uu2
	    + 2.0*g(I23,i)*uu2*uu3
	    + g(I33,i)*uu3*uu3;
	
	  const Real W = std::sqrt(1.0 + usq);
	  const Real alpha = std::sqrt(-1.0/gi(I00,i));

	  const Real u_t = -alpha*W + g(I01,i)*uu1 + g(I02,i)*uu2 + g(I03,i)*uu3;

	  const Real h = 1.0 + gamma_gas/(gamma_gas-1.0)*press/rho;

	  const Real rho_ut = pmb->phydro->u(IDN,k,j,i);

	  if (rho_ut <= 0.0) continue;
	
	  const Real dm = vol(i)*rho_ut;

	  // Three definitions of specific ejecta energy
	  const Real egeom = -u_t - 1.0;
	  const Real ebern = -h*u_t - 1.0;
	  const Real ebind = CellSpecificBindingEnergy(pmb, k,j,i, gi);
	  
	  // --------------------------------------------------
	  // Ejecta masses and energies
	  // --------------------------------------------------

	  if (egeom > 0.0) {
	    cache(MGEOM) += dm;
	    cache(EGEOM) += dm*egeom;
	  }
	
	  if (ebern > 0.0) {
	    cache(MBERN) += dm;
	    cache(EBERN) += dm*ebern;
	  }
	
	  if (ebind > 0.0) {
	    cache(MBIND) += dm;
	    cache(EBIND) += dm*ebind;
	  }
	
	  // --------------------------------------------------
	  // Homology diagnostics:
	  // geodesically-unbound outward matter
	  // --------------------------------------------------
	
	  if (egeom > 0.0) {
	  
	    // dr/dt; beta^r = 0 in the present metric
	    const Real vr = alpha*uu1/W;
	  
	    if (vr > 0.0) {
	      cache(MGEOM_OUT) += dm;

	      const Real v2 = 1.0 - 1.0/(W*W);

	      if (v2 > 0.0) {
		const Real chiP = press/(rho*v2);
		
		if (chiP < chiP_cut) {
		  cache(MBAL) += dm;
		}
	      }

	      const Real r = pmb->pcoord->x1v(i);

	      const Real texp = r/vr;

	      cache(MTEXP)  += dm*texp;
	      cache(MTEXP2) += dm*texp*texp;
	    
	      cache(MR2)  += dm*r*r;
	      cache(MRVR) += dm*r*vr;
	      cache(MVR2) += dm*vr*vr;

	      cache(RMAX) = std::max(cache(RMAX), r);
	    }
	  }

	}
      }
    }
    
    cache(HCYCLE) = static_cast<Real>(ncycle);
  }
  

  Real HistoryEjectaMassGeom(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);

    const Real M =
      pmb->ruser_meshblock_data[2](MGEOM);
    
    return M / pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEjectaEnergyGeom(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    const Real E =
      pmb->ruser_meshblock_data[2](EGEOM);
    
    return E * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      / pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryEjectaMassBern(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);

    const Real M =
      pmb->ruser_meshblock_data[2](MBERN);
    
    return M / pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEjectaEnergyBern(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    const Real E =
      pmb->ruser_meshblock_data[2](EBERN);
    
    return E * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      / pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryEjectaMassBind(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);

    const Real M =
      pmb->ruser_meshblock_data[2](MBIND);
    
    return M / pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEjectaEnergyBind(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    const Real E =
      pmb->ruser_meshblock_data[2](EBIND);
    
    return E * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      / pmb->pmy_mesh->punit->erg_code;

  }
  
  Real HistoryEjectaMaxRadius(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real r_max = pmb->ruser_meshblock_data[2](RMAX);
    return r_max*pmb->pmy_mesh->punit->code_length_cgs;
    
  }

  Real HistoryEjectaBallisticMass(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MBAL);
    return mass / pmb->pmy_mesh->punit->solar_mass_code;
    
  }

  Real HistoryEjectaMassExpTime(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MTEXP);
    return mass / pmb->pmy_mesh->punit->solar_mass_code * pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryEjectaMassExpTimeSq(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MTEXP2);
    return mass / pmb->pmy_mesh->punit->solar_mass_code * SQR(pmb->pmy_mesh->punit->code_time_cgs);
  }

  Real HistoryEjectaMassRadiusSq(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MR2);
    return mass / pmb->pmy_mesh->punit->solar_mass_code * SQR(pmb->pmy_mesh->punit->code_length_cgs);
  }

  Real HistoryEjectaMassRadiusVelocity(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MRVR);
    return mass / pmb->pmy_mesh->punit->solar_mass_code * pmb->pmy_mesh->punit->code_length_cgs / pmb->pmy_mesh->punit->speed_of_light_code;
    // Units of Msun cm c
  }

  Real HistoryEjectaMassGeomOut(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MGEOM_OUT);
    return mass / pmb->pmy_mesh->punit->solar_mass_code;
    // Units of Msun cm c
  }

  Real HistoryEjectaMassVelocitySq(MeshBlock *pmb, int iout) {
    UpdateEjectaHistoryCache(pmb);
    
    Real mass = pmb->ruser_meshblock_data[2](MVR2);
    return mass / pmb->pmy_mesh->punit->solar_mass_code / SQR(pmb->pmy_mesh->punit->speed_of_light_code);
    // Units of Msun c^2
  }


  bool EjectaReachedRadius(MeshBlock *pmb, Real rstop) {

    // This MeshBlock does not reach rstop.
    if (pmb->pcoord->x1v(pmb->ie) < rstop) {
      return false;
    }

    AthenaArray<Real> g, gi;
    g.NewAthenaArray(NMETRIC, pmb->ncells1);
    gi.NewAthenaArray(NMETRIC, pmb->ncells1);

    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {

	pmb->pcoord->CellMetric(
				k, j, pmb->is, pmb->ie, g, gi);

	// Search from outside inward.
	for (int i=pmb->ie; i>=pmb->is; --i) {

	  const Real r = pmb->pcoord->x1v(i);

	  if (r < rstop) {
	    break;
	  }

	  const Real uu1 = pmb->phydro->w(IVX,k,j,i);
	  const Real uu2 = pmb->phydro->w(IVY,k,j,i);
	  const Real uu3 = pmb->phydro->w(IVZ,k,j,i);

	  const Real usq =
            g(I11,i)*uu1*uu1
	    + 2.0*g(I12,i)*uu1*uu2
	    + 2.0*g(I13,i)*uu1*uu3
	    + g(I22,i)*uu2*uu2
	    + 2.0*g(I23,i)*uu2*uu3
	    + g(I33,i)*uu3*uu3;

	  const Real W = std::sqrt(1.0 + usq);
	  const Real alpha = std::sqrt(-1.0/gi(I00,i));

	  const Real u_t =
            -alpha*W
	    + g(I01,i)*uu1
	    + g(I02,i)*uu2
	    + g(I03,i)*uu3;

	  const Real egeom = -u_t - 1.0;
	  const Real vr = alpha*uu1/W;

	  if (egeom > 0.0 && vr > 0.0) {
	    return true;
	  }
	}
      }
    }

    return false;
  }

  Real HistoryOuterMassFlux(MeshBlock *pmb, int iout) {
    Real mdot = pmb->ruser_meshblock_data[3](MDOT_OUT);
    // make it units of Msun s^-1
    return mdot
      /pmb->pmy_mesh->punit->solar_mass_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryMdotInBC(MeshBlock *pmb, int iout) {
    
    Real mdot = pmb->ruser_meshblock_data[3](MDOT_IN_BC);
    return mdot/pmb->pmy_mesh->punit->solar_mass_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryMdotOutBC(MeshBlock *pmb, int iout) {

    Real mdot = pmb->ruser_meshblock_data[3](MDOT_OUT_BC);
    return mdot/pmb->pmy_mesh->punit->solar_mass_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryEdotInBC(MeshBlock *pmb, int iout) {
    
    Real edot = pmb->ruser_meshblock_data[3](EDOT_IN_BC);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryEdotOutBC(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](EDOT_OUT_BC);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryMdotInNum(MeshBlock *pmb, int iout) {
    
    Real mdot = pmb->ruser_meshblock_data[3](MDOT_IN_NUM);
    return mdot/pmb->pmy_mesh->punit->solar_mass_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryMdotOutNum(MeshBlock *pmb, int iout) {

    Real mdot = pmb->ruser_meshblock_data[3](MDOT_OUT_NUM);
    return mdot/pmb->pmy_mesh->punit->solar_mass_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryEdotInNum(MeshBlock *pmb, int iout) {
    
    Real edot = pmb->ruser_meshblock_data[3](EDOT_IN_NUM);
    return edot
      *SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryEdotOutNum(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](EDOT_OUT_NUM);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }
  
  Real HistoryMassOut(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_OUT);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryMassInBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_IN_BC);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryMassOutBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_OUT_BC);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEnergyInBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](ENERGY_IN_BC);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryEnergyOutBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](ENERGY_OUT_BC);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryMassInNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_IN_NUM);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryMassOutNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_OUT_NUM);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryEnergyInNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](ENERGY_IN_NUM);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryEnergyOutNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](ENERGY_OUT_NUM);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryBdotInBC(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](BDOT_IN_BC);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryBdotInNum(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](BDOT_IN_NUM);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }
  
  Real HistoryBdotOutBC(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](BDOT_OUT_BC);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryBdotOutNum(MeshBlock *pmb, int iout) {

    Real edot = pmb->ruser_meshblock_data[3](BDOT_OUT_NUM);
    return edot
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code
      /pmb->pmy_mesh->punit->code_time_cgs;
  }

  Real HistoryBindInBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](BIND_IN_BC);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryBindOutBC(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](BIND_OUT_BC);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryBindInNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](BIND_IN_NUM);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryBindOutNum(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real energy = pm->ruser_mesh_data[2](BIND_OUT_NUM);
    return energy
      * SQR(pmb->pmy_mesh->punit->speed_of_light_code)
      /pmb->pmy_mesh->punit->erg_code;
  }

  Real HistoryMassInjOut(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_INJ_OUT);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryMassInjIn(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_INJ_IN);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }

  Real HistoryRemovedMass(MeshBlock *pmb, int iout) {
    Mesh *pm = pmb->pmy_mesh;
    Real mass = pm->ruser_mesh_data[2](MASS_REMOVED);
    return mass/pmb->pmy_mesh->punit->solar_mass_code;
  }
  
} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // cell size
  if (mesh_size.x1min <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "Logarithmic x1 grid requires x1min > 0." << std::endl;
    ATHENA_ERROR(msg);
  }

  EnrollUserMeshGenerator(X1DIR, LogarithmicX1);

  const Real x1rat =
    std::pow(mesh_size.x1max/mesh_size.x1min,
	     1.0/static_cast<Real>(mesh_size.nx1));
  
  if (Globals::my_rank == 0) {
    std::cout << "effective x1rat = "
	      << std::setprecision(17) << x1rat << std::endl;
  }

  gamma_gas = pin->GetReal("hydro", "gamma");

  Real m_bh_code = 0.0;
  Real ang_bh_code = 0.0;
  
  if(!restart_){
    // Collapsed progenitor
    std::string filename = pin->GetString("problem", "progenitor_file");
    
    const bool initialize_with_mass =
      pin->GetOrAddBoolean("coord", "initialize_with_mass", true);
    
    ProgenitorProfile progenitor = ReadProgenitorProfile(filename);

    Real t0=0.0;
    if(initialize_with_mass){
      Real bh_mass_init = pin->GetReal("problem", "initial_bh_mass");
      t0 = GetTimeFromBlackHoleMass(progenitor, bh_mass_init);
      if(Globals::my_rank==0) std::cout << "BH mass is specified to " << bh_mass_init << " Msun." << std::endl;
    }else{
      t0 = pin->GetReal("problem", "t0");
      if(Globals::my_rank==0) std::cout << "Time is specified." << std::endl;
    }
    if(Globals::my_rank==0) std::cout << "t0 = " << t0 << " s." << std::endl;

    collapsed = CollapseProgenitor(progenitor, t0);
    
    Real rin_code = mesh_size.x1min;
    Real rin_cgs = rin_code*punit->code_length_cgs;
    Real M_inner_cgs = collapsed.EnclosedMass(rin_cgs);
    Real J_inner_cgs = 0.0;
    
    m_bh_code = Constants::grav_const_cgs * M_inner_cgs
      / SQR(Constants::speed_of_light_cgs)
      / punit->code_length_cgs;
    
    ang_bh_code = J_inner_cgs
      * std::pow(Constants::grav_const_cgs / SQR(Constants::speed_of_light_cgs)
		 / punit->code_length_cgs, 2);
    
    if (Globals::my_rank == 0) {
      printf("Black hole mass (cgs,code unit)=%15.7e, %15.7e\n",M_inner_cgs, m_bh_code);
    }
  }
  
  // pin->SetReal("coord", "m", m_bh_code);
  // pin->SetReal("coord", "j", ang_bh_code);
  // m_bh_code = pin->GetReal("coord", "m");

  // std::cout << std::setprecision(17)
  //           << "m_bh_code = " << m_bh_code
  //           << "  m_from_pin = " << pin->GetReal("coord","m")
  //           << std::endl;

  if (mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::user) {
    if(Globals::my_rank == 0) {
      std::cout << "Injecting BC with outflow_file..." << std::endl;
    }
    // injecting BC
    std::string fname = pin->GetString("problem", "outflow_file");
    poutflow = new OutflowBoundaryData(fname);

    const Real rin_cgs  = mesh_size.x1min * punit->code_length_cgs;
    const Real rext_cgs = poutflow->GetExtractionRadius();

    const Real relerr = std::abs(rin_cgs - rext_cgs) / rext_cgs;
    
    if (relerr > 1.0e-10) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
	  << "Inner boundary radius is inconsistent with the "
	  << "outflow extraction radius." << std::endl
	  << "Rin   = " << rin_cgs << " cm" << std::endl
	  << "r_ext = " << rext_cgs << " cm" << std::endl
	  << "relative difference = " << relerr << std::endl;
      ATHENA_ERROR(msg);
    }
    
    if(Globals::my_rank==0)poutflow->Analyze(gamma_gas);
    
    timescale_cut = pin->GetReal("problem", "timescale_cut");
    W_max = pin->GetReal("problem", "W_max");
    
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1,
			       InjectionInnerX1);
  }else{
    if(Globals::my_rank == 0) {
      std::cout << "ix1_bc = "
		<< GetBoundaryString(mesh_bcs[BoundaryFace::inner_x1])
		<< std::endl;
    }
  }

  // if (mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::user) {
  //   EnrollUserBoundaryFunction(BoundaryFace::outer_x1,
  // 			       DiodeOuterX1);
  // }
  // output
  AllocateUserHistoryOutput(48);

  int iout=0;
  EnrollUserHistoryOutput(iout, HistoryBlackHoleMass, "m_bh",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBlackHoleMassAccretionRate, "mdot_bh",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassOut, "mass_out",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryOuterMassFlux, "mdot_out",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassGeom, "ejecta_mass(geom)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaEnergyGeom, "ejecta_energy(geom)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassBern, "ejecta_mass(bern)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaEnergyBern, "ejecta_energy(bern)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassBind, "ejecta_mass(bind)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaEnergyBind, "ejecta_energy(bind)",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMaxRadius, "maximum_radius",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMInnerBlocks, "M_inner_block",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEInnerBlocks, "E_inner_block",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryJInnerBlocks, "J_inner_block",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMdotInBC, "mdot_in_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMdotOutBC, "mdot_out_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEdotInBC, "edot_in_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEdotOutBC, "edot_out_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMdotInNum, "mdot_in_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMdotOutNum, "mdot_out_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEdotInNum, "edot_in_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEdotOutNum, "edot_out_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassInBC, "mass_in_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassOutBC, "mass_out_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEnergyInBC, "energy_in_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEnergyOutBC, "energy_out_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassInNum, "mass_in_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassOutNum, "mass_out_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEnergyInNum, "energy_in_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEnergyOutNum, "energy_out_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBdotInBC, "bdot_in_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBdotOutBC, "bdot_out_BC",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBdotInNum, "bdot_in_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBdotOutNum, "bdot_out_num",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBindInBC, "bind_in_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBindOutBC, "bind_out_BC",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBindInNum, "bind_in_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryBindOutNum, "bind_out_num",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassInjIn, "mass_inj_in",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryMassInjOut, "mass_inj_out",
                          UserHistoryOperation::max);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassGeomOut, "mgeom_out",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaBallisticMass, "mejecta_ballistic",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassExpTime, "mej_texp",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassExpTimeSq, "mej_texp^2",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassRadiusSq, "mej_r^2",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassRadiusVelocity, "mej_r_vr",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryEjectaMassVelocitySq, "mej_vr^2",
                          UserHistoryOperation::sum);
  ++iout;
  EnrollUserHistoryOutput(iout, HistoryRemovedMass, "removed_mass",
                          UserHistoryOperation::max);


  //
  AllocateRealUserMeshDataField(3);

  // BH mass
  ruser_mesh_data[0].NewAthenaArray(1);
  // BH spin
  ruser_mesh_data[1].NewAthenaArray(1);

  ruser_mesh_data[0](0) = m_bh_code;
  ruser_mesh_data[1](0) = ang_bh_code;

  // Integrated boundary fluxes
  ruser_mesh_data[2].NewAthenaArray(NBOUNDARY_FLUX_INTEGRAL);
  ruser_mesh_data[2].ZeroClear();
}

void Mesh::UserWorkInLoop(void) {

  Real flux_sum[NBOUNDARY_FLUX_INTEGRAL] = {};

  for (int b = 0; b < nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);

    const AthenaArray<Real> &cache = pmb->ruser_meshblock_data[3];

    flux_sum[MASS_IN_BC]  += cache(MDOT_IN_BC);
    flux_sum[ENERGY_IN_BC]  += cache(EDOT_IN_BC);
    flux_sum[MASS_IN_NUM] += cache(MDOT_IN_NUM);
    flux_sum[ENERGY_IN_NUM] += cache(EDOT_IN_NUM);

    flux_sum[MASS_OUT_BC]  += cache(MDOT_OUT_BC);
    flux_sum[ENERGY_OUT_BC]  += cache(EDOT_OUT_BC);
    flux_sum[MASS_OUT_NUM] += cache(MDOT_OUT_NUM);
    flux_sum[ENERGY_OUT_NUM] += cache(EDOT_OUT_NUM);

    flux_sum[BIND_IN_BC] += cache(BDOT_IN_BC);
    flux_sum[BIND_IN_NUM] += cache(BDOT_IN_NUM);
    flux_sum[BIND_OUT_BC] += cache(BDOT_OUT_BC);
    flux_sum[BIND_OUT_NUM] += cache(BDOT_OUT_NUM);

    flux_sum[MASS_OUT] += cache(MDOT_OUT);

    flux_sum[MASS_INJ_IN] += cache(MDOT_INJ_IN);
    flux_sum[MASS_INJ_OUT] += cache(MDOT_INJ_OUT);
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
		flux_sum,
		NBOUNDARY_FLUX_INTEGRAL,
		MPI_ATHENA_REAL,
		MPI_SUM,
		MPI_COMM_WORLD);
#endif
  
  for (int n = 0; n < NBOUNDARY_FLUX_INTEGRAL; ++n) {
    ruser_mesh_data[2](n) += dt * flux_sum[n];
  }


  // write criteria for termination
  bool terminate = false;

  const Real rstop = 0.95 * mesh_size.x1max;
   
  int reached = 0;

  for (int b=0; b<nblocal; ++b) {
    if (EjectaReachedRadius(my_blocks(b), rstop)) {
      reached = 1;
      break;
    }
  }
  

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &reached, 1,
                MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  if (reached) terminate = true;
  
  if (terminate) {
    tlim = time + dt;

    if (Globals::my_rank == 0) {
      std::cout << "Termination condition satisfied." << std::endl;
    }
  }
  
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {

  const Real mass_threshold =
      pin->GetOrAddReal("problem",
                        "inner_block_mass_threshold", 0.01);

  Real M_inner = 0.0;
  Real E_inner = 0.0;
  Real r_ej_max = 0.0;
  int ejecta_reached = 0;

  const Real rstop = 0.95 * mesh_size.x1max;

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);

    M_inner += HistoryMInnerBlocks(pmb, 0);
    E_inner += HistoryEInnerBlocks(pmb, 0);

    r_ej_max = std::max(
        r_ej_max, HistoryEjectaMaxRadius(pmb, 0));

    if (EjectaReachedRadius(pmb, rstop)) {
      ejecta_reached = 1;
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &M_inner, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &E_inner, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &r_ej_max, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ejecta_reached, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  const int add_outer = ejecta_reached ? 1 : 0;
  
  // A first-pass mass-based recommendation.
  const int remove_inner = (nrbx1 > 1 && M_inner < mass_threshold) ? 1 : 0;

  if (Globals::my_rank == 0) {
    std::cout << "\n=== Restart recommendation ===\n"
              << "time             = " << time * punit->code_time_cgs << " s\n"
              << "cycle            = " << ncycle << "\n"
              << "M_inner_block    = " << M_inner << " Msun\n"
              << "E_inner_block    = " << E_inner << " erg\n"
              << "r_ej_max         = " << r_ej_max * punit->code_length_cgs
              << " cm\n"
              << "r_outer          = " << mesh_size.x1max * punit->code_length_cgs
              << " cm\n"
              << "ejecta_reached   = " << ejecta_reached << "\n\n"
              << "<restart_mesh>\n"
              << "remove_inner_blocks = " << remove_inner << "\n"
              << "add_outer_blocks    = " << add_outer << "\n"
              << "==============================\n";
  }
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // model-specific output fields
  auto *pdrv = pmy_mesh->pmetric_driver;
  const int n_output = pdrv->NumModelOutputVariables();
  
  AllocateUserOutputVariables(6+n_output);
  
  SetUserOutputVariableName(0, "gtt");
  SetUserOutputVariableName(1, "grr");
  SetUserOutputVariableName(2, "alpha-1");
  SetUserOutputVariableName(3, "Lorentz-1");
  SetUserOutputVariableName(4, "u_t+1");
  SetUserOutputVariableName(5, "enthalpy-1");

  for(int n=0; n<n_output; ++n){
    SetUserOutputVariableName(6+n, pdrv->ModelOutputVariableName(n));
  }

  //
  AllocateRealUserMeshBlockDataField(4);
  
  // [0] for Phi
  ruser_meshblock_data[0].NewAthenaArray(ncells1 + 1);
  // [1] for Mgrav
  ruser_meshblock_data[1].NewAthenaArray(ncells1 + 1);

  // [2] ejecta/history diagnostics cache
  ruser_meshblock_data[2].NewAthenaArray(NEJECTA_CACHE);

  // [3] boundary flux diagnostics cache
  ruser_meshblock_data[3].NewAthenaArray(NBOUNDARY_FLUX_CACHE);
  
  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();
  ruser_meshblock_data[2].ZeroClear();
  ruser_meshblock_data[3].ZeroClear();
  
  ruser_meshblock_data[2](HCYCLE) = -1.0;


  // Name pscalars
#if NSCALARS > 0
  pscalars->SetScalarName(IYE,    "Ye");
  pscalars->SetScalarName(IXINJ,  "Xinj");
  pscalars->SetScalarName(IXPHYS, "Xphys");

  pscalars->SetScalarName(IXNEUT, "Xn");
  pscalars->SetScalarName(IXH1,   "XH1");
  pscalars->SetScalarName(IXHE3,  "XHe3");
  pscalars->SetScalarName(IXHE4,  "XHe4");
  pscalars->SetScalarName(IXC12,  "XC12");
  pscalars->SetScalarName(IXN14,  "XN14");
  pscalars->SetScalarName(IXO16,  "XO16");
  pscalars->SetScalarName(IXNE20, "XNe20");
  pscalars->SetScalarName(IXMG24, "XMg24");
  pscalars->SetScalarName(IXSI28, "XSi28");
  pscalars->SetScalarName(IXS32,  "XS32");
  pscalars->SetScalarName(IXAR36, "XAr36");
  pscalars->SetScalarName(IXCA40, "XCa40");
  pscalars->SetScalarName(IXTI44, "XTi44");
  pscalars->SetScalarName(IXCR48, "XCr48");
  pscalars->SetScalarName(IXCR56, "XCr56");
  pscalars->SetScalarName(IXFE52, "XFe52");
  pscalars->SetScalarName(IXFE54, "XFe54");
  pscalars->SetScalarName(IXFE56, "XFe56");
  pscalars->SetScalarName(IXNI56, "XNi56");
#endif
  }


void MeshBlock::UserWorkInLoop(void) {
}

void MeshBlock::UserWorkAfterHydroFlux(int stage, Real t_stage) {
  if(stage !=2 ) return;

  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ncells1 + 1);
  gi.NewAthenaArray(NMETRIC, ncells1 + 1);

  Real mdot_in_bc = 0.0;
  Real edot_in_bc = 0.0;
  Real mdot_in_num = 0.0;
  Real edot_in_num = 0.0;

  Real mdot_out_bc = 0.0;
  Real edot_out_bc = 0.0;
  Real mdot_out_num = 0.0;
  Real edot_out_num = 0.0;

  Real bdot_in_bc = 0.0;
  Real bdot_in_num = 0.0;
  Real bdot_out_bc = 0.0;
  Real bdot_out_num = 0.0;

  if (pbval->block_bcs[BoundaryFace::inner_x1] != BoundaryFlag::block){

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {

	pmetric->Face1Metric(k, j, is, is, g, gi);
	const int i = is;
	const Real 
	  &g_00 = g(I00,i), &g_01 = g(I01,i), &g_02 = g(I02,i), &g_03 = g(I03,i),
	  &g_10 = g(I01,i), &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
	  &g_20 = g(I02,i), &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
	  &g_30 = g(I03,i), &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
	const Real
	  &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
	  &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
	  &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
	  &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);
	
	const Real alpha = 1.0/std::sqrt(-g00);
	
	const Real area = pcoord->GetFace1Area(k,j,i);

	// primitives in the first ghost cell
	const Real rho = phydro->w(IDN,k,j,i-1);
	const Real pgas= phydro->w(IPR,k,j,i-1);
	const Real uu1 = phydro->w(IVX,k,j,i-1);
	const Real uu2 = phydro->w(IVY,k,j,i-1);
	const Real uu3 = phydro->w(IVZ,k,j,i-1);
	const Real u_sq = g_11*SQR(uu1) + 2.0*g_12*uu1*uu2 + 2.0*g_13*uu1*uu3
                        + g_22*SQR(uu2) + 2.0*g_23*uu2*uu3
	                + g_33*SQR(uu3);
	const Real gamma = std::sqrt(1.0 + u_sq);
	const Real egas = peos->EgasFromRhoP(rho, pgas);
	const Real wgas = rho + egas + pgas;
	// u^a
	const Real ucon0 = gamma/alpha;
	const Real ucon1 = uu1 - alpha*gamma*g01;
	const Real J1 = rho*ucon1;
	const Real T10 = wgas*ucon1*ucon0 + pgas*g01;

	const Real F_D_bc = alpha*J1;
	const Real F_tau_bc = alpha*alpha*T10 - F_D_bc;

	if(F_D_bc>0.0){
	  mdot_in_bc  += area * F_D_bc;
	  edot_in_bc  += area * F_tau_bc;
	  bdot_in_bc  += area * (alpha*(F_tau_bc + F_D_bc) - F_D_bc);
	}else{
	  mdot_out_bc  -= area * F_D_bc;
	  edot_out_bc  -= area * F_tau_bc;
	  bdot_out_bc  -= area * (alpha*(F_tau_bc + F_D_bc) - F_D_bc);
	}
	
	const Real F_D_num = phydro->flux[X1DIR](IDN,k,j,is);
	const Real F_tau_num = phydro->flux[X1DIR](IEN,k,j,is);
	if(F_D_num>0.0){
	  mdot_in_num += area * F_D_num;
	  edot_in_num += area * F_tau_num;
	  bdot_in_num += area * (alpha*(F_tau_num + F_D_num) - F_D_num);
	}else{
	  mdot_out_num -= area * F_D_num;
	  edot_out_num -= area * F_tau_num;
	  bdot_out_num -= area * (alpha*(F_tau_num + F_D_num) - F_D_num);
	}
      }
    }
  }

  //std::cout << "mdot = " << mdot_out_num << " " << mdot_in_num << std::endl;

  // outer boundary
  Real mdot_out = 0.0;
  if (pbval->block_bcs[BoundaryFace::outer_x1] != BoundaryFlag::block){

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
	const int i = ie; //index of outermost active cell
	const Real area = pcoord->GetFace1Area(k,j,i+1);
	const Real F_D_num = phydro->flux[X1DIR](IDN,k,j,i+1);
	mdot_out += area * F_D_num;
      }
    }
  }  
  
  AthenaArray<Real> &cache = ruser_meshblock_data[3];

  cache(MDOT_IN_BC) = mdot_in_bc;
  cache(EDOT_IN_BC) = edot_in_bc;
  cache(MDOT_IN_NUM) = mdot_in_num;
  cache(EDOT_IN_NUM) = edot_in_num;
  cache(MDOT_OUT_BC) = mdot_out_bc;
  cache(EDOT_OUT_BC) = edot_out_bc;
  cache(MDOT_OUT_NUM) = mdot_out_num;
  cache(EDOT_OUT_NUM) = edot_out_num;

  cache(BDOT_IN_BC) = bdot_in_bc;
  cache(BDOT_IN_NUM) = bdot_in_num;
  cache(BDOT_OUT_BC) = bdot_out_bc;
  cache(BDOT_OUT_NUM) = bdot_out_num;

  cache(MDOT_OUT) = mdot_out;
}


void MeshBlock::UserWorkAfterScalarFlux(int stage, Real t_stage) {
  if(stage !=2 ) return;

  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ncells1 + 1);
  gi.NewAthenaArray(NMETRIC, ncells1 + 1);

  Real mdot_inj_in = 0.0;
  Real mdot_inj_out = 0.0;

  if (pbval->block_bcs[BoundaryFace::inner_x1] != BoundaryFlag::block){

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
	
	const Real area = pcoord->GetFace1Area(k,j,is);
	
	const Real F_Dinj = pscalars->s_flux[X1DIR](IXINJ,k,j,is);
	if(F_Dinj>0.0){
	  mdot_inj_in += area * F_Dinj;
	}else{
	  mdot_inj_out -= area * F_Dinj;
	}
      }
    }
  }
  
  AthenaArray<Real> &cache = ruser_meshblock_data[3];
  cache(MDOT_INJ_IN) = mdot_inj_in;
  cache(MDOT_INJ_OUT) = mdot_inj_out;
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


  // model-specific output
  auto *pdrv = pmy_mesh->pmetric_driver;
  const int n_output = pdrv->NumModelOutputVariables();

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
	for (int n=0; n<n_output; ++n) {
	  user_out_var(6+n,k,j,i) = pdrv->ModelOutputVariable(this,n,k,j,i);
	}
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
      std::strcmp(COORDINATE_SYSTEM, "gr_dynamic_valencia") != 0 && 
      std::strcmp(COORDINATE_SYSTEM, "gr_dynamic") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_collapsar.cpp" << std::endl
        << "This problem generator requires Schwarzschild/Minkowski/gr_dynamic/gr_dynamic_valencia coordinates."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  
  Real rho_atmos = pin->GetReal("problem", "rho_atmos");
  Real press_atmos = pin->GetReal("problem", "press_atmos");
  int ind_first = collapsed.FirstCellIndex();
  
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real r = pcoord->x1v(i);
        Real rad_cgs = r * pmy_mesh->punit->code_length_cgs;
        
        int ind = collapsed.CellIndexFromRadius(rad_cgs);
        
        Real rho_cgs, press_cgs, jrot_cgs, uu1, uu2, uu3;
	Real ye, xphys, xinj, xprog[ProgenitorSpecies::NPROG_SPECIES];
        uu2 = 0.0;
        uu3 = 0.0;
	xinj = 0.0;
	xphys = 0.0;
        if( ind>=collapsed.ncell ){
          // outside the star atmosphere
          rho_cgs = rho_atmos;
          press_cgs = press_atmos;
          uu1 = 0.0;
	  jrot_cgs = 0.0;
	  ye = collapsed.ye[collapsed.ncell-1];
	  for(int n=0; n<ProgenitorSpecies::NPROG_SPECIES; ++n){
	    xprog[n] = collapsed.x[n][collapsed.ncell-1];
	  }
        }else if (ind==collapsed.ncell-1){
          // between the last cell center and stellar surface
          rho_cgs   = collapsed.rho[ind];
          press_cgs = collapsed.press[ind];
          uu1 = collapsed.ur[ind] / Constants::speed_of_light_cgs;
	  jrot_cgs = collapsed.jrot[ind];
	  ye = collapsed.ye[ind];
	  xphys = 1.0;
	  for(int n=0; n<ProgenitorSpecies::NPROG_SPECIES; ++n){
	    xprog[n] = collapsed.x[n][ind];
	  }
        }else if (ind<0){
          rho_cgs = collapsed.rho[ind_first];
          press_cgs= collapsed.press[ind_first];
          uu1 = collapsed.ur[ind_first]/Constants::speed_of_light_cgs;
	  jrot_cgs = collapsed.jrot[ind_first];
          ye = collapsed.ye[ind_first];
	  xphys = 1.0;
	  for(int n=0; n<ProgenitorSpecies::NPROG_SPECIES; ++n){
	    xprog[n] = collapsed.x[n][ind_first];
	  }
        } else {
          Real xx1 = (rad_cgs-collapsed.radius[ind])/(collapsed.radius[ind+1]-collapsed.radius[ind]);
          Real xx0 = 1.0-xx1;
          rho_cgs  = xx0*collapsed.rho[ind] + xx1*collapsed.rho[ind+1];
          press_cgs= xx0*collapsed.press[ind] + xx1*collapsed.press[ind+1];
          uu1 = (xx0*collapsed.ur[ind] + xx1*collapsed.ur[ind+1])/Constants::speed_of_light_cgs;
	  jrot_cgs = xx0*collapsed.jrot[ind] + xx1*collapsed.jrot[ind+1];
	  ye = xx0*collapsed.ye[ind] + xx1*collapsed.ye[ind+1];
	  xphys = 1.0;
	  for(int n=0; n<ProgenitorSpecies::NPROG_SPECIES; ++n){
	    xprog[n] = xx0*collapsed.x[n][ind] + xx1*collapsed.x[n][ind+1];
	  }
        }
        Real rho_code = rho_cgs/pmy_mesh->punit->code_density_cgs;
        Real press_code = press_cgs/pmy_mesh->punit->code_pressure_cgs;
	Real jrot_code = jrot_cgs /(pmy_mesh->punit->code_velocity_cgs*pmy_mesh->punit->code_length_cgs);
	// jrot is mass-shell-averaged value. jrot_local is (3/2)jrot*sin(theta)^2
	// Therefore, u^3 = g^33 u_3 = (3/2)*jrot/(r^2 h).
	Real h = 1.0+gamma_gas/(gamma_gas-1.0) * press_code/rho_code;
	uu3 = 1.5*jrot_code/(r*r*h);

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
        //printf("i, rho, press, uu1 = %5d %12.4e %12.4e %12.4e\n", i, rho_code, press_code, uu1);

#if NSCALARS > 0
	// pscalars
	pscalars->r(IYE,k,j,i) = ye;
	pscalars->r(IXINJ,k,j,i) = 0.0;
	pscalars->r(IXPHYS,k,j,i) = xphys;
	for (int n = IXNEUT; n < NSCALAR_REQUIRED; ++n) {
	  pscalars->r(n,k,j,i) = 0.0;
	}
	for (int n = 0; n < ProgenitorSpecies::NPROG_SPECIES; ++n) {
	  pscalars->r(prog_to_scalar[n],k,j,i) += xprog[n];
	}
#endif
      }
    }
  }

#if DYNAMIC_METRIC_ENABLED
  const bool initial_metric_update =
    pin->GetOrAddBoolean("coord", "initial_metric_update", true);
  
  if (!initial_metric_update) {
    std::stringstream msg;
    msg << "### FATAL ERROR in gr_collapsar::ProblemGenerator" << std::endl
	<< "Dynamic metric fresh start requires initial_metric_update=true."
	<< std::endl;
    ATHENA_ERROR(msg);
  }
#endif


#if !DYNAMIC_METRIC_ENABLED
  // Convert primitive -> conserved
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  bb.ZeroClear();

  peos->PrimitiveToConserved(
      phydro->w, bb, phydro->u, pcoord,
      is, ie, js, je, ks, ke);

#if NSCALARS > 0
  peos->PassiveScalarPrimitiveToConserved(
    pscalars->r, phydro->u, pscalars->s, pcoord,
    is, ie, js, je, ks, ke);
#endif // NSCALARS > 0
#endif // !DYNAMIC_METRIC_ENABLED  
}


void MeshBlock::InitializeAtmosphere(ParameterInput *pin) {
  const int il = is - NGHOST;
  const int iu = ie + NGHOST;
  const int jl = (block_size.nx2 > 1) ? js-NGHOST : js;
  const int ju = (block_size.nx2 > 1) ? je+NGHOST : je;
  const int kl = (block_size.nx3 > 1) ? ks-NGHOST : ks;
  const int ku = (block_size.nx3 > 1) ? ke+NGHOST : ke;

  const Real rho_cgs = pin->GetReal("problem", "rho_atmos");
  const Real press_cgs = pin->GetReal("problem", "press_atmos");

  const Real rho = rho_cgs / pmy_mesh->punit->code_density_cgs;
  const Real press = press_cgs / pmy_mesh->punit->code_pressure_cgs;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        phydro->w(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = press;
        phydro->w(IVX,k,j,i) = 0.0;
        phydro->w(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = 0.0;

        phydro->w1(IDN,k,j,i) = rho;
        phydro->w1(IPR,k,j,i) = press;
        phydro->w1(IVX,k,j,i) = 0.0;
        phydro->w1(IVY,k,j,i) = 0.0;
        phydro->w1(IVZ,k,j,i) = 0.0;

#if NSCALARS > 0
	pscalars->r(IYE,k,j,i)   = 0.5;
	pscalars->r(IXINJ,k,j,i) = 0.0;
	pscalars->r(IXPHYS,k,j,i) = 0.0;
	
	for (int n=IXNEUT; n<NSCALAR_REQUIRED; ++n) {
	  pscalars->r(n,k,j,i) = 0.0;
	}
	pscalars->r(IXHE4,k,j,i) = 1.0;
#endif

      }
    }
  }

#if !DYNAMIC_METRIC_ENABLED
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  bb.ZeroClear();

  peos->PrimitiveToConserved(
      phydro->w, bb, phydro->u, pcoord,
      is, ie, js, je, ks, ke);

#if NSCALARS > 0
  peos->PassiveScalarPrimitiveToConserved(
    pscalars->r, phydro->u, pscalars->s, pcoord,
    is, ie, js, je, ks, ke);
#endif // NSCALARS > 0
#endif // !DYNAMIC_METRIC_ENABLED  
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

  const Real G = Constants::grav_const_cgs;
  const Real c = Constants::speed_of_light_cgs;
  const Real Msun = Constants::solar_mass_cgs;
  
  const Real m_bh_cgs = bh_mass * Msun;

  // find radius for m=bh_mass  
  const int nface = progenitor.nface;
  int ilo, ihi;
  ilo = 0;
  ihi = nface-1;
  
  if(m_bh_cgs>=progenitor.mass_face[ihi]){
    std::stringstream msg;
    msg << "### FATAL ERROR in GetTimeFromBlackHoleMass" << std::endl
        << "Specified BH mass is larger than total stellar mass."
        << "BH mass = " << m_bh_cgs << "Mtot = " << progenitor.mass_face[ihi]
        << std::endl;
    ATHENA_ERROR(msg);
  }
  
  while(ihi-ilo>=2){
    int i=(ihi+ilo)/2;
    if ( progenitor.mass_face[i] <= m_bh_cgs ){
      ilo = i;
    }else{
      ihi = i;
    }
  }
  // Radius of enclosed mass = m_bh_cgs
  Real w1 = (std::log(m_bh_cgs)-std::log(progenitor.mass_face[ilo]))
    /(std::log(progenitor.mass_face[ihi])-std::log(progenitor.mass_face[ilo]));
  Real w0 = 1.0 - w1;
  Real r_m0 = std::exp( w0*std::log(progenitor.radius_face[ilo]) + w1*std::log(progenitor.radius_face[ihi]) );

  // derive collapse parameter eta for which r(m) = 2*m.
  Real eta = std::acos(4.0*G*m_bh_cgs/(c*c*r_m0) - 1.0);

  // derive sound crossing time
  Real t_m0 = 0.0;
  for(int i=0; i<ilo; ++i){
    Real dr = progenitor.radius_face[i+1] - progenitor.radius_face[i];
    Real dt = dr/progenitor.csound[i];
    t_m0 += dt;
  }
  t_m0 += (r_m0 - progenitor.radius_face[ilo])/ progenitor.csound[ilo];

  Real t_ff = std::sqrt(r_m0*r_m0*r_m0/(8.0*G*m_bh_cgs)) * (eta + std::sin(eta));
  Real tau = t_m0 + t_ff;

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

  collapsed.ye = progenitor.ye;
  collapsed.x = progenitor.x;

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

  const Real tmax = pmb->pmy_mesh->poutflow->GetTimeMax();
  const Real tmin = pmb->pmy_mesh->poutflow->GetTimeMin();
  
  const Real vv_max = 1.0 - 1.0/(W_max*W_max);
  const Real time_cgs = time * code_time_cgs;

  const Real time_interp = std::max(tmin, std::min(time_cgs, tmax));
  const Real s = std::max(0.0, std::min((time_cgs - tmax)/timescale_cut, 1.0));

  const Real w_inj = 1.0 - s*s*(3.0 - 2.0*s);
  
  // Real fac = 1.0;
  // Real time_interp = time_cgs;

  // // Avoid extrapolation
  // // should determine how we will deal with after running out of data.
  // if      ( time_cgs > pmb->pmy_mesh->poutflow->GetTimeMax() ){
  //   time_interp = pmb->pmy_mesh->poutflow->GetTimeMax();
  //   fac = std::exp(- (time_cgs-pmb->pmy_mesh->poutflow->GetTimeMax())/timescale_cut );
  // }else if( time_cgs < pmb->pmy_mesh->poutflow->GetTimeMin() ){
  //   time_interp = pmb->pmy_mesh->poutflow->GetTimeMin();
  //   fac = 1.0;
  // }
  
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
      for (int ng=1; ng<=ngh; ++ng) {
	
        const int i = il - ng;
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
        
        prim(IDN,k,j,i) = state.rho/code_rho_cgs;
        prim(IPR,k,j,i) = state.press/code_press_cgs;
        prim(IVX,k,j,i) = uu1;
        prim(IVY,k,j,i) = uu2;
        prim(IVZ,k,j,i) = uu3;

	for (int n=0; n<NHYDRO; ++n) {
	  Real diode = prim(n,k,j,il);
	  if (n == IVX) {
	    diode = std::min(diode, 0.0);
	  }
	  prim(n,k,j,i) = w_inj*prim(n,k,j,i) + (1.0-w_inj)*diode;
	}
	
	
#if NSCALARS > 0
	pmb->pscalars->r(IYE,k,j,i)   = state.ye;
	pmb->pscalars->r(IXINJ,k,j,i) = 1.0;
	pmb->pscalars->r(IXPHYS,k,j,i) = 1.0;
	for (int n=IXNEUT; n<NSCALAR_REQUIRED; ++n) {
	  pmb->pscalars->r(n,k,j,i) = 0.0;
	}
#endif
      }
    }
  }
}


// void DiodeInnerX1(MeshBlock *pmb, Coordinates *pco,
//                   AthenaArray<Real> &prim, FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju,
//                   int kl, int ku, int ngh) {

//   for (int k=kl; k<=ku; ++k) {
//     for (int j=jl; j<=ju; ++j) {

//       const Real ur =
//           std::min(prim(IVX,k,j,il), 0.0);

//       for (int i=il-ngh; i<il; ++i) {

//         // Copy all primitive hydro variables
//         for (int n=0; n<NHYDRO; ++n) {
//           prim(n,k,j,i) = prim(n,k,j,il);
//         }

//         // Do not allow outward radial velocity
//         prim(IVX,k,j,i) = ur;

// #if NSCALARS > 0
//         // Copy passive scalar mass fractions
//         for (int n=0; n<NSCALARS; ++n) {
//           pmb->pscalars->r(n,k,j,i) =
//               pmb->pscalars->r(n,k,j,il);
//         }
// #endif
//       }
//     }
//   }
// }

// void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco,
//                   AthenaArray<Real> &prim, FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju,
//                   int kl, int ku, int ngh) {

//   for (int k=kl; k<=ku; ++k) {
//     for (int j=jl; j<=ju; ++j) {

//       const Real ur = std::max(prim(IVX,k,j,iu), 0.0);
      
//       for (int i=iu+1; i<=iu+ngh; ++i) {

//         // Copy all primitive hydro variables
//         for (int n=0; n<NHYDRO; ++n) {
//           prim(n,k,j,i) = prim(n,k,j,iu);
//         }

//         // Do not allow inward radial velocity
//         prim(IVX,k,j,i) = ur;

// #if NSCALARS > 0
//         // Copy passive scalar mass fractions
//         for (int n=0; n<NSCALARS; ++n) {
//           pmb->pscalars->r(n,k,j,i) = pmb->pscalars->r(n,k,j,iu);
//         }
// #endif
//       }
//     }
//   }
// }



Real LogarithmicX1(Real x, RegionSize rs) {
  return rs.x1min
    * std::pow(rs.x1max/rs.x1min, x);
}

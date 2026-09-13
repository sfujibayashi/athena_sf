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

namespace HelmholtzConstants {
  const Real ssol=5.6704e-5, amu=1.66053878283e-24, h=6.6260689633e-27;
  const Real qe=4.8032042712e-10, avo=6.0221417930e23, clight=2.99792458e10,
                    kerg=1.380650424e-16;
  const Real asol=4.0*ssol/clight, light2=clight*clight;
  const Real asoli3=asol/3.0;
  const Real kergavo=kerg*avo, sioncon = (2.0 * PI * amu * kerg)/(h*h);

  //const Real qe=4.8032042712e-10;
  const Real esqu=qe*qe;

  const Real eV_to_erg = 1.602176634e-12, MeV_to_erg = 1.0e6 * eV_to_erg;
} // namespace HelmholtzConstants

Real threshold;
Real gamma_gas;

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
  }
  return;
}

namespace helm{

  //HelmTable* phelm = nullptr;
  int i_ye = -1;
  int i_ytot = -1;
  int i_temp = -1;
  int i_mexc = -1;
}

namespace tracer {
  constexpr int i_r0 = 4;
}

namespace uov{
  int i_entr = 0;
  int i_temp = 1;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(2);
  
  //if (!helm::phelm) helm::phelm = new HelmTable(pin, ptable);
  helm::i_ye   = pin->GetInteger("hydro", "helm_ye_index");
  helm::i_ytot = pin->GetInteger("hydro", "helm_ytot_index");
  helm::i_temp = pin->GetInteger("hydro", "helm_temp_index");
  helm::i_mexc = pin->GetInteger("hydro", "helm_mexc_index");
  std::cout << helm::i_ye <<" "<<helm::i_ytot<<" "<<helm::i_temp<<" "<<helm::i_mexc<<std::endl;

  if (tracer::i_r0 >= NSCALARS) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::InitUserMeshBlockData\n"
        << "Not enough passive scalars for initial-radius tracer.\n";
    ATHENA_ERROR(msg);
  }

  if (tracer::i_r0 == helm::i_ye   ||
      tracer::i_r0 == helm::i_ytot ||
      tracer::i_r0 == helm::i_temp ||
      tracer::i_r0 == helm::i_mexc) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::InitUserMeshBlockData\n"
        << "Initial-radius tracer index is equal to either of i_ye, i_ytot, i_temp, i_mexc, which are booked in Helmholtz EOS.\n";
    ATHENA_ERROR(msg);
  }
  
  gamma_gas = pin->GetReal("hydro", "gamma");
  
}

void MeshBlock::UserWorkInLoop(void) {
  
  Real v1_max = 0.0;
  int i_vmax = -1;
  Real rho_vmax = 0.0;
  Real mom_vmax = 0.0;
  Real cs2_max = 0.0;
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
	Real rho = phydro->u(IDN,k,j,i);
	Real temp= pscalars->r(helm::i_temp,k,j,i);
	Real ye  = pscalars->r(helm::i_ye,k,j,i);
	Real ytot= pscalars->r(helm::i_ytot,k,j,i);
        Real abar= 1.0/ytot;
#if HELMHOLTZ_EOS_ENABLED
	AthenaArray<Real> out;
	out.NewAthenaArray(8);
	peos->HelmLookupRhoT(rho, temp, ye, abar, out);
        Real asq = out(4);
#else
        Real asq = gamma_gas*phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);;
#endif
        if (cs2_max < asq){
          cs2_max=asq;
        }
        Real mom = phydro->u(IM1,k,j,i);
        Real v1 = mom/rho;
        if (v1_max < std::abs(v1)){
          v1_max = std::abs(v1);
          i_vmax = i;
          rho_vmax = rho;
          mom_vmax = mom;
        }
      }
    }
  }
  printf("v1max=%12.4e i=%d rho=%12.4e mom=%12.4e cs2max=%12.4e\n",
         v1_max, i_vmax, rho_vmax, mom_vmax, cs2_max);

  bool isok;
  isok = true;
  // printf("%12.4e %12.4e\n", phydro->u(IEN,0,0,0), phydro->u(IEN,0,0,0)/phydro->u(IDN,0,0,0));
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
	for (int n=0; n<NHYDRO; n++) {
	  if(not isfinite(phydro->u(n,k,j,i))){
	    //printf("hydro, %3d%3d\n",n,i);
	    isok=false;
	  }
	}
	
	for (int n=0; n<NSCALARS; n++) {
	  if(not isfinite(pscalars->r(n,k,j,i))){
	    //printf("scalar, %3d%3d\n",n,i);
	    isok=false;
	  }
	}
	
      }
    }
  }

  // int k=ks;
  // int j=js;
  // for(int i=is; i<=ie; i++) {
  //   printf("%4d %12.4e",i,pcoord->x1v(i));
  //   for (int n=0; n<NWAVE; n++) {
  //     printf("%12.4e",phydro->u(n,j,k,i));
  //   }
  //   for (int n=0; n<NSCALARS; n++) {
  //     printf("%12.4e",pscalars->s(n,j,k,i));
  //   }

  //   printf("  ");
  //   for (int n=0; n<NWAVE; n++) {
  //     printf("%12.4e",phydro->w(n,j,k,i));
  //   }
  //   for (int n=0; n<NSCALARS; n++) {
  //     printf("%12.4e",pscalars->r(n,j,k,i));
  //   }
  //   printf("\n");
  // }


  if(not isok){
    std::stringstream msg;
    msg << "NaN detected.\n";
    ATHENA_ERROR(msg);
  }

}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      for(int i=is; i<=ie; i++) {
	Real rho = phydro->u(IDN,k,j,i);
	Real temp= pscalars->r(helm::i_temp,k,j,i);
	Real ye  = pscalars->r(helm::i_ye,k,j,i);
	Real ytot= pscalars->r(helm::i_ytot,k,j,i);
        Real abar= 1.0/ytot;
#if HELMHOLTZ_EOS_ENABLED
	AthenaArray<Real> out;
	out.NewAthenaArray(8);
	peos->HelmLookupRhoT(rho, temp, ye, abar, out);
	Real entr = out(7);
#else
        Real entr = phydro->w(IPR,k,j,i)/std::pow(phydro->w(IDN,k,j,i),gamma_gas);
#endif
	//printf("%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",rho,temp,ye,abar,entr,out(0),out(2),out(5));
	user_out_var(uov::i_entr,k,j,i) = entr;
	user_out_var(uov::i_temp,k,j,i) = temp;
      }
    }
  }

}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rout = pin->GetReal("problem", "rout");
  Real rin  = pin->GetReal("problem", "rin");
  
  Real Mej  = pin->GetReal("problem", "ejecta_mass");
  Real Vinf = pin->GetReal("problem", "asymptotic_velocity");
  Real mexc = pin->GetReal("problem", "mass_excess");
  Real abar = pin->GetReal("problem", "mass_number");
  Real ye   = pin->GetReal("problem", "ye");

  Real ratio_atmos   = pin->GetReal("problem", "ratio_atmos");
  
  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast_helmholtz.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  Real vol = 4.0*PI/3.0 * (pow(rout,3.0) - pow(rin,3.0));
  Real Eint = 0.5*Mej*Vinf*Vinf;
  
  Real rho_ej = Mej / vol;
  Real eth_ej = Eint/ vol;
  
  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        Real rho, eth, v1;
        if (rin < rad and rad < rout) {
	  rho = rho_ej;
	  eth = eth_ej;
	  v1  = (rad/rout)*Vinf;
	} else {
	  rho = rho_ej*ratio_atmos;
	  eth = eth_ej*ratio_atmos;
	  v1  = 0.0;
        }
	//v1 = 0.0;
	
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho*v1;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = eth + 0.5*rho*v1*v1;

      }
    }
  }
  
  
  if (NSCALARS > 0) {

    // first initialize all pscalars
    for (int n=0; n<NSCALARS; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pscalars->s(n,k,j,i) = 0.0;
          }
        }
      }
    }

    // Ye
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  pscalars->s(helm::i_ye,k,j,i) = ye * phydro->u(IDN,k,j,i);
	}
      }
    }

    // Ytot = 1/Abar
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
          Real ytot = 1.0 / abar;
	  pscalars->s(helm::i_ytot,k,j,i) = ytot * phydro->u(IDN,k,j,i);
	}
      }
    }

    // mass excess per baryon (MeV)
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  pscalars->s(helm::i_mexc,k,j,i) = mexc * phydro->u(IDN,k,j,i);
	}
      }
    }

    // initial radius as a tracer
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  Real rad;	    
	  if        (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
	    rad = sqrt( SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)) );
	  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0){
	    rad = sqrt( SQR(pcoord->x1v(i)) + SQR(pcoord->x3v(k)) );
	  } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
	    rad = pcoord->x1v(i);
	  }
	  pscalars->s(tracer::i_r0,k,j,i) = rad*phydro->u(IDN,k,j,i);
	}
      }
    }
  }
  
#if MASS_EXCESS_ENERGY_ENABLED
  {
    using namespace HelmholtzConstants;
    // add mass-excess contribution
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  phydro->u(IEN,k,j,i) += pscalars->s(helm::i_mexc,k,j,i)*MeV_to_erg*avo;
	}
      }
    }
  }
#endif

  {
    using namespace HelmholtzConstants;
    
    // derive temperature
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  Real rho  = phydro->u(IDN,k,j,i);
          Real ke = 0.5 / rho *
            (  + SQR(phydro->u(IM1,k,j,i))
               + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)) );

	  Real egas = phydro->u(IEN,k,j,i) - ke;
	  Real s_cell[NSCALARS];
	  for (int n=0; n<NSCALARS; ++n) {
	    s_cell[n] = pscalars->s(n,k,j,i);
	  }
	  s_cell[helm::i_temp] = 1e10*phydro->u(IDN,k,j,i);
	  // std::cout << i << " " << rho << " "<< egas  << std::endl;
	  // std::cout << s[0]/rho << " " << s[1]/rho << " "<< s[2]/rho <<" "<< s[3]/rho << std::endl;
	  // std::cout << pscalars->s(0,k,j,i)/rho << " " << pscalars->s(1,k,j,i)/rho << " "<< pscalars->s(2,k,j,i)/rho <<" "<< pscalars->s(3,k,j,i)/rho << std::endl;
#if HELMHOLTZ_EOS_ENABLED	  
	  Real temp = peos->TempFromRhoEg(rho, egas, s_cell);
#else
          Real temp = 0.0;
#endif
	  pscalars->s(helm::i_temp,k,j,i) = temp * phydro->u(IDN,k,j,i);
	  printf("%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n", pcoord->x1v(i), rho, temp, egas, egas/rho, s_cell[helm::i_mexc]*MeV_to_erg*avo/rho);
	  
	}
      }
    }
  }


}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  MeshBlock *pmb = my_blocks(0);

  // analysis - check shape of the spherical blast wave
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pmb->phydro->w, 4, IPR, 1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast_helmholtz.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic=is; ic<=ie; ic++)
    if (pmb->pcoord->x1f(ic) > x1_0) break;
  ic--;
  for (jc=pmb->js; jc<=pmb->je; jc++)
    if (pmb->pcoord->x2f(jc) > x2_0) break;
  jc--;
  for (kc=pmb->ks; kc<=pmb->ke; kc++)
    if (pmb->pcoord->x3f(kc) > x3_0) break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o=0; o<=6; o++) {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1) ios=-10;
    else if (o == 2) ios =  10;
    else if (o == 3) jos = -10;
    else if (o == 4) jos =  10;
    else if (o == 5) kos = -10;
    else if (o == 6) kos =  10;
    for (int d=0; d<6; d++) {
      Real pmax = 0.0;
      int imax(0), jmax(0), kmax(0);
      if (d == 0) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i>=is; i--) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 1) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i<=ie; i++) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 2) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j>=js; j--) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 3) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j<=je; j++) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 4) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k>=ks; k--) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      } else { // if (d == 5) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k<=ke; k++) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      Real x1m = pmb->pcoord->x1v(imax);
      Real x2m = pmb->pcoord->x2v(jmax);
      Real x3m = pmb->pcoord->x3v(kmax);
      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else {  // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if (rad > rmax) rmax = rad;
      if (rad < rmin) rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  Real  x1c = pmb->pcoord->x1v(ic);
  Real dx1c = pmb->pcoord->dx1f(ic);
  Real  x2c = pmb->pcoord->x2v(jc);
  Real dx2c = pmb->pcoord->dx2f(jc);
  Real dx3c = pmb->pcoord->dx3f(kc);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
    }
    std::fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    std::fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    std::fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    std::fclose(pfile);
  }
  return;
}


// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  if (pmb->pmy_mesh->f3) {
    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
          maxeps = std::max(maxeps, eps);
        }
      }
    }
  } else if (pmb->pmy_mesh->f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
        maxeps = std::max(maxeps, eps);
      }
    }
  } else {
    return 0;
  }

  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1;
  return 0;
}

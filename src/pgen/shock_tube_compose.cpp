//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file shock_tube_compose.cpp
//! \brief 1D shock tube test for the CompOSE tabulated EOS.
//========================================================================================

#include <cmath>
#include <sstream>

#include "../athena.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/compose_table.hpp"


#if NSCALARS < 1
#error "shock_tube_compose requires at least one passive scalar for Ye."
#endif


namespace {

  constexpr Real amu_cgs = 1.66053906660e-24;
  constexpr Real mevfm3_to_cgs = 1.602176634e33;


  //----------------------------------------------------------------------------------------
  // Return pressure in code units from rho [code], T [MeV], and Ye.
  
  Real PressureFromRhoTYe(EquationOfState *peos,
                          Real rho, Real temp, Real ye,
                          Real rho_unit, Real egas_unit) {
    Real rho_cgs = rho * rho_unit;
    Real nb = rho_cgs / amu_cgs * 1.0e-39;
    
    Real log_nb = std::log10(nb);
    Real log_t = std::log10(temp);
    
    Real log_p =
      peos->ptable->table3d.interpolate(ECLOGP, log_nb, ye, log_t);
    
    return std::pow(10.0, log_p) * mevfm3_to_cgs / egas_unit;
  }

  //----------------------------------------------------------------------------------------
  // Return temperature [MeV] from rho [code], P [code], and Ye.

  Real TemperatureFromRhoPYe(EquationOfState *peos,
                             Real rho, Real pres, Real ye,
                             Real rho_unit, Real egas_unit) {
    Real rho_cgs = rho * rho_unit;
    Real nb = rho_cgs / amu_cgs * 1.0e-39;

    Real log_nb = std::log10(nb);
    Real log_p  = std::log10(pres * egas_unit / mevfm3_to_cgs);

    Real log_t_lo, log_t_hi;
    peos->ptable->table3d.GetX1lim(log_t_lo, log_t_hi);

    Real f_lo =
      peos->ptable->table3d.interpolate(ECLOGP, log_nb, ye, log_t_lo) - log_p;
    Real f_hi =
      peos->ptable->table3d.interpolate(ECLOGP, log_nb, ye, log_t_hi) - log_p;

    if (f_lo * f_hi > 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in TemperatureFromRhoPYe" << std::endl
          << "Pressure is outside the temperature range of the EOS table."
          << std::endl
          << "rho = " << rho
          << ", pres = " << pres
          << ", Ye = " << ye << std::endl;
      ATHENA_ERROR(msg);
    }

    for (int n = 0; n < 80; ++n) {
      Real log_t_mid = 0.5 * (log_t_lo + log_t_hi);

      Real f_mid =
        peos->ptable->table3d.interpolate(
                                          ECLOGP, log_nb, ye, log_t_mid) - log_p;

      if (std::abs(f_mid) < 1.0e-12 ||
          std::abs(log_t_hi - log_t_lo) < 1.0e-12) {
        return std::pow(10.0, log_t_mid);
      }

      if ((f_lo <= 0.0 && f_mid >= 0.0) ||
          (f_lo >= 0.0 && f_mid <= 0.0)) {
        log_t_hi = log_t_mid;
      } else {
        log_t_lo = log_t_mid;
        f_lo = f_mid;
      }
    }

    return std::pow(10.0, 0.5 * (log_t_lo + log_t_hi));
  }
  
} // namespace

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//! \brief Allocate user-defined output variable for temperature.
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(4);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "entr");
  SetUserOutputVariableName(2, "Yn");
  SetUserOutputVariableName(3, "Yp");
}


//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//! \brief Store temperature [MeV] in user-defined output variable.
//========================================================================================

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real rho_unit =
      pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0);
  Real egas_unit =
      pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0);

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real rho  = phydro->w(IDN,k,j,i);
        Real pres = phydro->w(IPR,k,j,i);
        Real ye   = pscalars->r(0,k,j,i);

        Real temp =
          TemperatureFromRhoPYe(peos, rho, pres, ye,
                                rho_unit, egas_unit);

        Real rho_cgs = rho * rho_unit;
        Real nb = rho_cgs / amu_cgs * 1.0e-39;
        
        Real log_nb = std::log10(nb);
        Real log_t  = std::log10(temp);

        user_out_var(0,k,j,i) = temp;
        user_out_var(1,k,j,i) =
          peos->ptable->table3d.interpolate(ECENT, log_nb, ye, log_t);
        user_out_var(2,k,j,i) =
          peos->ptable->table3d.interpolate(ECYNEUT, log_nb, ye, log_t);
        user_out_var(3,k,j,i) =
          peos->ptable->table3d.interpolate(ECYPROT, log_nb, ye, log_t);
        
      }
    }
  }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initialize a 1D shock tube using the CompOSE EOS.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real xshock = pin->GetOrAddReal("problem", "xshock", 0.0);

  if (xshock < pmy_mesh->mesh_size.x1min ||
      xshock > pmy_mesh->mesh_size.x1max) {
    std::stringstream msg;
    msg << "### FATAL ERROR in shock_tube_compose" << std::endl
        << "xshock lies outside the x1 domain." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Left state
  Real rho_l  = pin->GetReal("problem", "rho_l");
  Real temp_l = pin->GetReal("problem", "temp_l");
  Real ye_l   = pin->GetReal("problem", "ye_l");
  Real vx_l   = pin->GetOrAddReal("problem", "vx_l", 0.0);

  // Right state
  Real rho_r  = pin->GetReal("problem", "rho_r");
  Real temp_r = pin->GetReal("problem", "temp_r");
  Real ye_r   = pin->GetReal("problem", "ye_r");
  Real vx_r   = pin->GetOrAddReal("problem", "vx_r", 0.0);

  Real rho_unit =
      pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0);
  Real egas_unit =
      pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0);

  // Pressure from rho, T, Ye.
  Real pres_l =
      PressureFromRhoTYe(peos, rho_l, temp_l, ye_l,
                         rho_unit, egas_unit);

  Real pres_r =
      PressureFromRhoTYe(peos, rho_r, temp_r, ye_r,
                         rho_unit, egas_unit);

  // Primitive scalar arrays supplied to the EOS.
  Real r_l[(NSCALARS > 0) ? NSCALARS : 1] = {};
  Real r_r[(NSCALARS > 0) ? NSCALARS : 1] = {};

  r_l[0] = ye_l;
  r_r[0] = ye_r;

  // Internal energy density in code units.
  Real egas_l =
      peos->EgasFromRhoP(rho_l, pres_l, r_l);

  Real egas_r =
      peos->EgasFromRhoP(rho_r, pres_r, r_r);

  // Initialize conserved variables.
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real rho, vx, ye, egas;

        if (pcoord->x1v(i) < xshock) {
          rho  = rho_l;
          vx   = vx_l;
          ye   = ye_l;
          egas = egas_l;
        } else {
          rho  = rho_r;
          vx   = vx_r;
          ye   = ye_r;
          egas = egas_r;
        }

        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * vx;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        phydro->u(IEN,k,j,i) =
            egas + 0.5 * rho * vx * vx;

        // Conserved scalar s = rho * Ye.
        pscalars->s(0,k,j,i) = rho * ye;
      }
    }
  }
}

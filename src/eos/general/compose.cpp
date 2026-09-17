//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file compose.cpp
//! \brief CompOSE tabulated EOS for the general EOS framework.
//========================================================================================

#include <cmath>
#include <limits>
#include <sstream>

#include "../../athena.hpp"
#include "../../parameter_input.hpp"
#include "../../utils/interp_table.hpp"
#include "../eos.hpp"
#include "../../utils/compose_table.hpp"

namespace {

  constexpr Real amu_cgs = 1.66053906660e-24;
  constexpr Real amu_mev = 931.49410242;
  constexpr Real clight = 2.99792458e10;
  constexpr Real mevfm3_to_cgs = 1.602176634e33;
  
  Real fixed_ye = -1.0;
  int i_ye = -1;
  
  
  // Convert rho [g cm^-3] to baryon number density [fm^-3].
  // The Athena++ rest-mass convention is rho = amu * nb.
  inline Real RhoToNb(Real rho) {
    return rho / amu_cgs * 1.0e-39;
  }
  

  // Return Ye from conserved scalar variables.
  // Conserved scalars are s_n = rho * r_n.
  inline Real YeFromConserved(Real rho, const Real *s) {
    if (i_ye >= 0) {
      if (s == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in YeFromConserved" << std::endl
            << "Composition scalars are required for CompOSE EOS."
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return s[i_ye] / rho;
    }
    return fixed_ye;
  }


  // Return Ye from primitive scalar variables.
  inline Real YeFromPrimitive(const Real *r) {
    if (i_ye >= 0) {
      if (r == nullptr) {
        std::stringstream msg;
        msg << "### FATAL ERROR in YeFromPrimitive" << std::endl
            << "Composition scalars are required for CompOSE EOS."
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return r[i_ye];
    }
    return fixed_ye;
  }


  // Check nb and Ye against the tabulated domain.
  inline void CheckTableCoordinates(EosTable *ptable, Real log_nb, Real ye) {
    Real log_nb_min, log_nb_max;
    Real ye_min, ye_max;

    ptable->table3d.GetX3lim(log_nb_min, log_nb_max);
    ptable->table3d.GetX2lim(ye_min, ye_max);

    if (!std::isfinite(log_nb) ||
        log_nb < log_nb_min || log_nb > log_nb_max) {
      std::stringstream msg;
      msg << "### FATAL ERROR in CompOSE EOS" << std::endl
          << "Baryon density is outside table range." << std::endl
          << "log10(nb) = " << log_nb
          << ", range = [" << log_nb_min << ", " << log_nb_max << "]"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    if (!std::isfinite(ye) || ye < ye_min || ye > ye_max) {
      std::stringstream msg;
      msg << "### FATAL ERROR in CompOSE EOS" << std::endl
          << "Ye is outside table range." << std::endl
          << "Ye = " << ye
          << ", range = [" << ye_min << ", " << ye_max << "]"
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }


  // Invert one table quantity along the temperature direction.
  //
  // table_var = 0 : log10(P)
  // table_var = 1 : log10(e_tot)
  //
  // target has the same (logarithmic) representation as the table quantity.
  Real FindLogTemperature(EosTable *ptable, int table_var,
                          Real log_nb, Real ye, Real target) {
    Real log_t_lo, log_t_hi;
    ptable->table3d.GetX1lim(log_t_lo, log_t_hi);

    Real f_lo =
      ptable->table3d.interpolate(table_var, log_nb, ye, log_t_lo) - target;
    Real f_hi =
      ptable->table3d.interpolate(table_var, log_nb, ye, log_t_hi) - target;

    if (f_lo == 0.0) return log_t_lo;
    if (f_hi == 0.0) return log_t_hi;

    if (f_lo * f_hi > 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in FindLogTemperature" << std::endl
          << "Target quantity is outside the temperature range of the EOS table."
          << std::endl
          << "target = " << target << std::endl
          << "table values at Tmin/Tmax = "
          << f_lo + target << ", " << f_hi + target << std::endl;
      ATHENA_ERROR(msg);
    }

    for (int n = 0; n < 80; ++n) {
      Real log_t_mid = 0.5 * (log_t_lo + log_t_hi);
      Real f_mid =
        ptable->table3d.interpolate(table_var, log_nb, ye, log_t_mid) - target;

      if (std::abs(f_mid) < 1.0e-12 ||
          std::abs(log_t_hi - log_t_lo) < 1.0e-12) {
        return log_t_mid;
      }

      if ((f_lo <= 0.0 && f_mid >= 0.0) ||
          (f_lo >= 0.0 && f_mid <= 0.0)) {
        log_t_hi = log_t_mid;
        f_hi = f_mid;
      } else {
        log_t_lo = log_t_mid;
        f_lo = f_mid;
      }
    }

    return 0.5 * (log_t_lo + log_t_hi);
  }

} // namespace

//----------------------------------------------------------------------------------------
//! \brief Return gas pressure from rho, internal energy density, and composition.

Real EquationOfState::PresFromRhoEg(Real rho, Real egas, Real *s) {
  Real ye = YeFromConserved(rho, s);

  // code density -> g cm^-3 -> fm^-3
  Real rho_cgs = rho * rho_unit_;
  Real nb = RhoToNb(rho_cgs);
  Real log_nb = std::log10(nb);

  CheckTableCoordinates(ptable, log_nb, ye);

  // Athena++ egas excludes the rest-mass contribution rho c^2.
  // Convert it to MeV fm^-3, then reconstruct total energy density
  // using the convention rho = amu * nb.
  Real eint =
      egas * egas_unit_ / mevfm3_to_cgs;
  Real etot = eint + amu_mev * nb;

  if (!std::isfinite(etot) || etot <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::PresFromRhoEg" << std::endl
        << "Non-positive total energy density." << std::endl
        << "rho = " << rho << ", egas = " << egas
        << ", etot = " << etot << std::endl;
    ATHENA_ERROR(msg);
  }

  Real log_etot = std::log10(etot);

  Real log_t = FindLogTemperature(ptable, ECLOGE, log_nb, ye, log_etot);

  Real log_p = ptable->table3d.interpolate(ECLOGP, log_nb, ye, log_t);
  Real pres_cgs = std::pow(10.0, log_p) * mevfm3_to_cgs;

  return pres_cgs * inv_egas_unit_;
}


//----------------------------------------------------------------------------------------
//! \brief Return internal energy density from rho, pressure, and composition.

Real EquationOfState::EgasFromRhoP(Real rho, Real pres, Real *r) {
  Real ye = YeFromPrimitive(r);

  Real rho_cgs = rho * rho_unit_;
  Real nb = RhoToNb(rho_cgs);
  Real log_nb = std::log10(nb);

  CheckTableCoordinates(ptable, log_nb, ye);

  Real pres_mevfm3 =
      pres * egas_unit_ / mevfm3_to_cgs;

  if (!std::isfinite(pres_mevfm3) || pres_mevfm3 <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::EgasFromRhoP" << std::endl
        << "Non-positive pressure." << std::endl
        << "rho = " << rho << ", pres = " << pres << std::endl;
    ATHENA_ERROR(msg);
  }

  Real log_p = std::log10(pres_mevfm3);

  Real log_t = FindLogTemperature(ptable, ECLOGP, log_nb, ye, log_p);

  Real log_etot = ptable->table3d.interpolate(ECLOGE, log_nb, ye, log_t);
  Real etot = std::pow(10.0, log_etot);

  // Remove the amu rest-mass contribution to obtain the Athena++ gas energy.
  Real eint = etot - amu_mev * nb;

  return eint * mevfm3_to_cgs * inv_egas_unit_;
}


//----------------------------------------------------------------------------------------
//! \brief Return adiabatic sound speed squared.

Real EquationOfState::AsqFromRhoP(Real rho, Real pres, const Real *r) {
  Real ye = YeFromPrimitive(r);

  Real rho_cgs = rho * rho_unit_;
  Real nb = RhoToNb(rho_cgs);
  Real log_nb = std::log10(nb);

  CheckTableCoordinates(ptable, log_nb, ye);

  Real pres_mevfm3 =
      pres * egas_unit_ / mevfm3_to_cgs;

  if (!std::isfinite(pres_mevfm3) || pres_mevfm3 <= 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::AsqFromRhoP" << std::endl
        << "Non-positive pressure." << std::endl
        << "rho = " << rho << ", pres = " << pres << std::endl;
    ATHENA_ERROR(msg);
  }

  Real log_p = std::log10(pres_mevfm3);

  Real log_t = FindLogTemperature(ptable, ECLOGP, log_nb, ye, log_p);

  // table variable 2 is cs/c
  Real cs = ptable->table3d.interpolate(ECCS, log_nb, ye, log_t);

  Real cs2_cgs = cs * cs * clight * clight;

  return cs2_cgs * inv_vsqr_unit_;
}


//----------------------------------------------------------------------------------------
//! \brief Initialize CompOSE EOS parameters.

void EquationOfState::InitEosConstants(ParameterInput *pin) {
  if (pin->DoesParameterExist("hydro", "compose_ye")) {
    fixed_ye = pin->GetReal("hydro", "compose_ye");
  }

  if (pin->DoesParameterExist("hydro", "compose_ye_index")) {
    i_ye = pin->GetInteger("hydro", "compose_ye_index");

    if (i_ye < 0 || i_ye >= NSCALARS) {
      std::stringstream msg;
      msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
          << "hydro/compose_ye_index must be between 0 and NSCALARS-1 ("
          << NSCALARS << ")." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  if (fixed_ye < 0.0 && i_ye < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState::InitEosConstants" << std::endl
        << "Either hydro/compose_ye or hydro/compose_ye_index must be specified."
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

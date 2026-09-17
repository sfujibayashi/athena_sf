//========================================================================================
//! \file compose_test.cpp
//! \brief Test problem for the CompOSE tabulated EOS.
//========================================================================================

#include <cmath>
#include <iostream>

#include "../athena.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/compose_table.hpp"


namespace {

constexpr Real amu_cgs = 1.66053906660e-24;
constexpr Real amu_mev = 931.49410242;
constexpr Real clight = 2.99792458e10;
constexpr Real mevfm3_to_cgs = 1.602176634e33;

} // namespace

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = my_blocks(0);
  EquationOfState *peos = pmb->peos;

  Real rho0 = pin->GetReal("problem", "rho");
  Real temp0= pin->GetReal("problem", "temp");
  Real ye0  = pin->GetReal("problem", "ye");

  Real rho_unit = pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0);
  Real egas_unit = pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0);
  Real vsqr_unit = egas_unit / rho_unit;

  // Code density -> physical density -> nb [fm^-3]
  Real rho_cgs = rho0 * rho_unit;
  Real nb = rho_cgs / amu_cgs * 1.0e-39;

  Real log_nb = std::log10(nb);
  Real log_t = std::log10(temp0);

  Real r_cell[(NSCALARS > 0) ? NSCALARS : 1] = {};
  Real s_cell[(NSCALARS > 0) ? NSCALARS : 1] = {};

  r_cell[0] = ye0;
  s_cell[0] = rho0 * ye0;

  // Direct lookup from the loaded table
  Real log_p =
    peos->ptable->table3d.interpolate(ECLOGP, log_nb, ye0, log_t);
  Real pres0 =
    std::pow(10.0, log_p) * mevfm3_to_cgs / egas_unit;
  Real egas0 = 
    peos->EgasFromRhoP(rho0, pres0, r_cell);
  
  Real max_rho_err = 0.0;
  Real max_ye_err  = 0.0;
  Real max_e_err = 0.0;

  for (int b = 0; b < nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real rho = pmb->phydro->u(IDN,k,j,i);
          Real ye  = pmb->pscalars->s(0,k,j,i) / rho;
          Real etot = pmb->phydro->u(IEN,k,j,i);
          
          max_rho_err = std::max(
              max_rho_err, std::abs((rho - rho0) / rho0));

          max_ye_err = std::max(
              max_ye_err, std::abs((ye - ye0) / ye0));

          max_e_err = std::max(
              max_e_err, std::abs((etot - egas0) / egas0));
        }
      }
    }
  }

  std::cout << "max relative rho error = "
            << max_rho_err << std::endl;
  std::cout << "max relative Ye error  = "
            << max_ye_err << std::endl;
  std::cout << "max relative Egas error  = "
            << max_e_err << std::endl;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rho = pin->GetReal("problem", "rho");
  Real temp = pin->GetReal("problem", "temp");
  Real ye = pin->GetReal("problem", "ye");

  Real rho_unit = pin->GetOrAddReal("hydro", "eos_rho_unit", 1.0);
  Real egas_unit = pin->GetOrAddReal("hydro", "eos_egas_unit", 1.0);
  Real vsqr_unit = egas_unit / rho_unit;

  // Code density -> physical density -> nb [fm^-3]
  Real rho_cgs = rho * rho_unit;
  Real nb = rho_cgs / amu_cgs * 1.0e-39;

  Real log_nb = std::log10(nb);
  Real log_t = std::log10(temp);

  // Direct lookup from the loaded table
  Real log_p =
      peos->ptable->table3d.interpolate(ECLOGP, log_nb, ye, log_t);
  Real log_etot =
      peos->ptable->table3d.interpolate(ECLOGE, log_nb, ye, log_t);
  Real cs =
      peos->ptable->table3d.interpolate(ECCS, log_nb, ye, log_t);

  Real pres_table =
      std::pow(10.0, log_p) * mevfm3_to_cgs / egas_unit;

  Real etot =
      std::pow(10.0, log_etot);

  Real egas_table =
      (etot - amu_mev * nb) * mevfm3_to_cgs / egas_unit;

  Real asq_table =
      cs * cs * clight * clight / vsqr_unit;

  // Primitive and conserved scalar representations
  Real r_cell[(NSCALARS > 0) ? NSCALARS : 1] = {};
  Real s_cell[(NSCALARS > 0) ? NSCALARS : 1] = {};

  r_cell[0] = ye;
  s_cell[0] = rho * ye;

  // Test EOS interface
  Real egas_eos =
      peos->EgasFromRhoP(rho, pres_table, r_cell);

  Real pres_eos =
      peos->PresFromRhoEg(rho, egas_table, s_cell);

  Real asq_eos =
      peos->AsqFromRhoP(rho, pres_table, r_cell);

  std::cout << "CompOSE EOS test" << std::endl;
  std::cout << "rho  = " << rho_cgs << " g/cm^3" << std::endl;
  std::cout << "nb   = " << nb << " fm^-3" << std::endl;
  std::cout << "T    = " << temp << " MeV" << std::endl;
  std::cout << "Ye   = " << ye << std::endl;

  std::cout << std::endl;
  std::cout << "P(table)   = " << pres_table << std::endl;
  std::cout << "P(EOS)     = " << pres_eos << std::endl;
  std::cout << "rel. error = "
            << (pres_eos - pres_table) / pres_table << std::endl;

  std::cout << std::endl;
  std::cout << "e(table)   = " << egas_table << std::endl;
  std::cout << "e(EOS)     = " << egas_eos << std::endl;
  std::cout << "rel. error = "
            << (egas_eos - egas_table) / egas_table << std::endl;

  std::cout << std::endl;
  std::cout << "asq(table) = " << asq_table << std::endl;
  std::cout << "asq(EOS)   = " << asq_eos << std::endl;
  std::cout << "rel. error = "
            << (asq_eos - asq_table) / asq_table << std::endl;


  // Initialize a uniform stationary state
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = egas_eos;

#if NSCALARS > 0
        pscalars->s(0,k,j,i) = rho * ye;
#endif
      }
    }
  }
}

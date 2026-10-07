//========================================================================================
//! \file valencia_flux_transform_test.cpp
//! \brief Unit test for Coordinates::FluxToGlobal1/2/3 in gr_dynamic_valencia.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../metric/metric.hpp"
#include "../parameter_input.hpp"

#if !GENERAL_RELATIVITY
#error "valencia_flux_transform_test requires GENERAL_RELATIVITY"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "valencia_flux_transform_test is hydro-only"
#endif

namespace {

struct LocalFlux {
  Real Jx, Txt, Txx, Txy, Txz;
};

struct MetricState {
  Real mass, psi, dm;
  const char *name;
};

inline Real ScaledError(Real x, Real xref) {
  return std::abs(x-xref) /
         std::max(static_cast<Real>(1.0), std::abs(xref));
}

inline void Check(const char *dir, const char *name, int k, int j, int i,
                  Real got, Real ref, Real tol,
                  Real &max_err, bool &failed, std::stringstream &log) {
  const Real err = ScaledError(got, ref);
  max_err = std::max(max_err, err);
  if (err > tol) {
    failed = true;
    log << std::setprecision(16)
        << "  FAIL " << dir << " " << name
        << " (k,j,i)=(" << k << "," << j << "," << i << ")"
        << " got=" << got << " ref=" << ref
        << " scaled_err=" << err << "\n";
  }
}

inline void MetricFactors(Real r, Real mass, Real psi, Real dm,
                          Real &alpha, Real &X) {
  const Real f   = 1.0 - 2.0*mass/r;
  const Real gtt = -f + 2.0*dm/r + 2.0*f*psi;
  const Real grr = 1.0/f + 2.0*dm/(r*f*f);

  if (!(f > 0.0) || !(gtt < 0.0) || !(grr > 0.0)) {
    std::stringstream msg;
    msg << "Invalid test metric at r=" << r
        << ": f=" << f << " gtt=" << gtt << " grr=" << grr;
    ATHENA_ERROR(msg);
  }
  alpha = std::sqrt(-gtt);
  X = std::sqrt(grr);
}

void SetMetricState(MeshBlock *pmb, const MetricState &s) {
  pmb->pmy_mesh->ruser_mesh_data[0](0) = s.mass;
  pmb->pmy_mesh->ruser_mesh_data[1](0) = 0.0;

  AthenaArray<Real> &psi = pmb->ruser_meshblock_data[0];
  AthenaArray<Real> &dm  = pmb->ruser_meshblock_data[1];
  for (int i=0; i<=pmb->ncells1; ++i) {
    psi(i) = s.psi;
    dm(i)  = s.dm;
  }
}

void SetX1(AthenaArray<Real> &f, int k, int j, int i, const LocalFlux &q) {
  f(IDN,k,j,i)=q.Jx;  f(IEN,k,j,i)=q.Txt;
  f(IM1,k,j,i)=q.Txx; f(IM2,k,j,i)=q.Txy; f(IM3,k,j,i)=q.Txz;
}

void SetX2(AthenaArray<Real> &f, int k, int j, int i, const LocalFlux &q) {
  // local (x,y,z) = (theta,phi,r)
  f(IDN,k,j,i)=q.Jx;  f(IEN,k,j,i)=q.Txt;
  f(IM2,k,j,i)=q.Txx; f(IM3,k,j,i)=q.Txy; f(IM1,k,j,i)=q.Txz;
}

void SetX3(AthenaArray<Real> &f, int k, int j, int i, const LocalFlux &q) {
  // local (x,y,z) = (phi,r,theta)
  f(IDN,k,j,i)=q.Jx;  f(IEN,k,j,i)=q.Txt;
  f(IM3,k,j,i)=q.Txx; f(IM1,k,j,i)=q.Txy; f(IM2,k,j,i)=q.Txz;
}

void RunCase(MeshBlock *pmb, const MetricState &state,
             const LocalFlux &lf, Real tol,
             Real &max_err, bool &failed, std::stringstream &log) {
  Coordinates *pco = pmb->pcoord;
  SetMetricState(pmb, state);

  AthenaArray<Real> flux, empty;
  flux.NewAthenaArray(NHYDRO, pmb->ncells3+1,
                      pmb->ncells2+1, pmb->ncells1+1);
  flux.ZeroClear();

  // x1: local x = r
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie+1; ++i) SetX1(flux,k,j,i,lf);

      pco->FluxToGlobal1(k,j,pmb->is,pmb->ie+1,
                         empty,empty,flux,empty,empty);

      for (int i=pmb->is; i<=pmb->ie+1; ++i) {
        const Real r = pco->x1f(i);
        const Real th = pco->x2v(j);
        Real alpha, X;
        MetricFactors(r,state.mass,state.psi,state.dm,alpha,X);

        Check("X1","D",k,j,i,flux(IDN,k,j,i),
              alpha/X*lf.Jx,tol,max_err,failed,log);
        Check("X1","tau",k,j,i,flux(IEN,k,j,i),
              alpha/X*(lf.Txt-lf.Jx),tol,max_err,failed,log);
        Check("X1","S1",k,j,i,flux(IM1,k,j,i),
              alpha*lf.Txx,tol,max_err,failed,log);
        Check("X1","S2",k,j,i,flux(IM2,k,j,i),
              alpha*r/X*lf.Txy,tol,max_err,failed,log);
        Check("X1","S3",k,j,i,flux(IM3,k,j,i),
              alpha*r*std::sin(th)/X*lf.Txz,tol,max_err,failed,log);

        Check("X1","q",k,j,i,
              pmb->pmetric->Face1SpatialDensitizationFactor(k,j,i),
              X,tol,max_err,failed,log);
      }
    }
  }

  // x2: local (x,y,z) = (theta,phi,r)
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je+1; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) SetX2(flux,k,j,i,lf);

      pco->FluxToGlobal2(k,j,pmb->is,pmb->ie,
                         empty,empty,flux,empty,empty);

      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real r = pco->x1v(i);
        const Real th = pco->x2f(j);
        const Real st = std::sin(th);
        Real alpha, X;
        MetricFactors(r,state.mass,state.psi,state.dm,alpha,X);

        Check("X2","D",k,j,i,flux(IDN,k,j,i),
              alpha/r*lf.Jx,tol,max_err,failed,log);
        Check("X2","tau",k,j,i,flux(IEN,k,j,i),
              alpha/r*(lf.Txt-lf.Jx),tol,max_err,failed,log);
        Check("X2","S1",k,j,i,flux(IM1,k,j,i),
              alpha*X/r*lf.Txz,tol,max_err,failed,log);
        Check("X2","S2",k,j,i,flux(IM2,k,j,i),
              alpha*lf.Txx,tol,max_err,failed,log);
        Check("X2","S3",k,j,i,flux(IM3,k,j,i),
              alpha*st*lf.Txy,tol,max_err,failed,log);

        Check("X2","q",k,j,i,
              pmb->pmetric->Face2SpatialDensitizationFactor(k,j,i),
              X,tol,max_err,failed,log);
      }
    }
  }

  // x3: local (x,y,z) = (phi,r,theta)
  for (int k=pmb->ks; k<=pmb->ke+1; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) SetX3(flux,k,j,i,lf);

      pco->FluxToGlobal3(k,j,pmb->is,pmb->ie,
                         empty,empty,flux,empty,empty);

      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real r = pco->x1v(i);
        const Real th = pco->x2v(j);
        const Real st = std::sin(th);
        Real alpha, X;
        MetricFactors(r,state.mass,state.psi,state.dm,alpha,X);

        Check("X3","D",k,j,i,flux(IDN,k,j,i),
              alpha/(r*st)*lf.Jx,tol,max_err,failed,log);
        Check("X3","tau",k,j,i,flux(IEN,k,j,i),
              alpha/(r*st)*(lf.Txt-lf.Jx),tol,max_err,failed,log);
        Check("X3","S1",k,j,i,flux(IM1,k,j,i),
              alpha*X/(r*st)*lf.Txy,tol,max_err,failed,log);
        Check("X3","S2",k,j,i,flux(IM2,k,j,i),
              alpha/st*lf.Txz,tol,max_err,failed,log);
        Check("X3","S3",k,j,i,flux(IM3,k,j,i),
              alpha*lf.Txx,tol,max_err,failed,log);

        Check("X3","q",k,j,i,
              pmb->pmetric->Face3SpatialDensitizationFactor(k,j,i),
              X,tol,max_err,failed,log);
      }
    }
  }

  std::cout << "  case: " << state.name << " done" << std::endl;
}

} // namespace

// MonopoleGravityDriver reads these before MeshBlocks are constructed.
void Mesh::InitUserMeshData(ParameterInput *pin) {
  AllocateRealUserMeshDataField(2);
  ruser_mesh_data[0].NewAthenaArray(1);  // BH mass
  ruser_mesh_data[1].NewAthenaArray(1);  // BH spin
  ruser_mesh_data[0](0) = 0.0;
  ruser_mesh_data[1](0) = 0.0;
}

// MonopoleGravityDriver uses these as Psi_face1 and delta_m_face1.
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(ncells1+1);
  ruser_meshblock_data[1].NewAthenaArray(ncells1+1);
  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real tol = pin->GetOrAddReal("problem","tol",1.0e-12);

  const LocalFlux lf = {
    static_cast<Real>(0.37),
    static_cast<Real>(1.11),
    static_cast<Real>(0.83),
    static_cast<Real>(-0.29),
    static_cast<Real>(0.47)
  };

  if (gid == 0) {
    const MetricState flat = {0.0,0.0,0.0,"flat spherical"};
    const MetricState nontrivial = {
      pin->GetOrAddReal("problem","test_mass",0.10),
      pin->GetOrAddReal("problem","test_psi",0.02),
      pin->GetOrAddReal("problem","test_dm",0.03),
      "nontrivial diagonal metric"
    };

    Real max_err = 0.0;
    bool failed = false;
    std::stringstream log;

    RunCase(this,flat,lf,tol,max_err,failed,log);
    RunCase(this,nontrivial,lf,tol,max_err,failed,log);

    std::cout << std::setprecision(16)
              << "Valencia FluxToGlobal unit test: max scaled error = "
              << max_err << std::endl;

    if (failed) {
      std::stringstream msg;
      msg << "### FATAL ERROR: Valencia FluxToGlobal unit test failed\n"
          << log.str();
      ATHENA_ERROR(msg);
    }

    std::cout << "Valencia FluxToGlobal unit test: PASS" << std::endl;
    SetMetricState(this,flat);
  }

  // Benign state for normal Athena++ initialization after the test.
  AthenaArray<Real> bb_cc;
  bb_cc.NewAthenaArray(3,ncells3,ncells2,ncells1);
  bb_cc.ZeroClear();

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->w(IDN,k,j,i)=1.0;
        phydro->w(IPR,k,j,i)=0.1;
        phydro->w(IVX,k,j,i)=0.0;
        phydro->w(IVY,k,j,i)=0.0;
        phydro->w(IVZ,k,j,i)=0.0;

        phydro->w1(IDN,k,j,i)=phydro->w(IDN,k,j,i);
        phydro->w1(IPR,k,j,i)=phydro->w(IPR,k,j,i);
        phydro->w1(IVX,k,j,i)=phydro->w(IVX,k,j,i);
        phydro->w1(IVY,k,j,i)=phydro->w(IVY,k,j,i);
        phydro->w1(IVZ,k,j,i)=phydro->w(IVZ,k,j,i);
      }
    }
  }

  peos->PrimitiveToConserved(phydro->w,bb_cc,phydro->u,pcoord,
                             is,ie,js,je,ks,ke);
}

void MeshBlock::InitializeAtmosphere(ParameterInput *pin) {
  return;
}

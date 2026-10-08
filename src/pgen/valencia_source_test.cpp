//========================================================================================
//! \file valencia_source_test.cpp
//! \brief Unit test for Valencia geometric source terms in gr_dynamic_valencia.
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
#include "../parameter_input.hpp"

#if !GENERAL_RELATIVITY
#error "valencia_source_test requires GENERAL_RELATIVITY"
#endif
#if MAGNETIC_FIELDS_ENABLED
#error "valencia_source_test is hydro-only"
#endif

namespace {

inline Real ScaledError(Real x, Real xref) {
  return std::abs(x-xref)/std::max(static_cast<Real>(1.0),std::abs(xref));
}

inline void Check(const char *cname, const char *vname,
                  int k, int j, int i, Real got, Real ref, Real tol,
                  Real &max_err, bool &failed, std::stringstream &log) {
  const Real err=ScaledError(got,ref);
  max_err=std::max(max_err,err);
  if (err>tol) {
    failed=true;
    log << std::setprecision(16)
        << "  FAIL " << cname << " " << vname
        << " (k,j,i)=(" << k << "," << j << "," << i << ")"
        << " got=" << got << " ref=" << ref
        << " scaled_err=" << err << "\n";
  }
}

void SetMetricState(MeshBlock *pmb, Real mass, Real psi, Real dm) {
  pmb->pmy_mesh->ruser_mesh_data[0](0)=mass;
  pmb->pmy_mesh->ruser_mesh_data[1](0)=0.0;
  AthenaArray<Real> &psi_f=pmb->ruser_meshblock_data[0];
  AthenaArray<Real> &dm_f =pmb->ruser_meshblock_data[1];
  for (int i=0; i<=pmb->ncells1; ++i) {
    psi_f(i)=psi;
    dm_f(i)=dm;
  }
}

void SetPrimitive(MeshBlock *pmb, Real rho, Real p,
                  Real u1, Real u2, Real u3) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        pmb->phydro->w(IDN,k,j,i)=rho;
        pmb->phydro->w(IPR,k,j,i)=p;
        pmb->phydro->w(IVX,k,j,i)=u1;
        pmb->phydro->w(IVY,k,j,i)=u2;
        pmb->phydro->w(IVZ,k,j,i)=u3;
        pmb->phydro->w1(IDN,k,j,i)=rho;
        pmb->phydro->w1(IPR,k,j,i)=p;
        pmb->phydro->w1(IVX,k,j,i)=u1;
        pmb->phydro->w1(IVY,k,j,i)=u2;
        pmb->phydro->w1(IVZ,k,j,i)=u3;
      }
    }
  }
}

void CopyActive(MeshBlock *pmb, const AthenaArray<Real> &src,
                AthenaArray<Real> &dst) {
  for (int n=0; n<NHYDRO; ++n)
    for (int k=pmb->ks; k<=pmb->ke; ++k)
      for (int j=pmb->js; j<=pmb->je; ++j)
        for (int i=pmb->is; i<=pmb->ie; ++i)
          dst(n,k,j,i)=src(n,k,j,i);
}

void P2C(MeshBlock *pmb, AthenaArray<Real> &bb_cc) {
  pmb->peos->PrimitiveToConserved(pmb->phydro->w,bb_cc,pmb->phydro->u,
                                  pmb->pcoord,
                                  pmb->is,pmb->ie,pmb->js,pmb->je,
                                  pmb->ks,pmb->ke);
}

void RunFlatStatic(MeshBlock *pmb, AthenaArray<Real> &bb_cc,
                   Real dt, Real tol, Real &max_err,
                   bool &failed, std::stringstream &log) {
  const Real rho=1.2;
  const Real p=0.17;
  SetMetricState(pmb,0.0,0.0,0.0);
  SetPrimitive(pmb,rho,p,0.0,0.0,0.0);
  P2C(pmb,bb_cc);

  AthenaArray<Real> u0,u1;
  u0.NewAthenaArray(NHYDRO,pmb->ncells3,pmb->ncells2,pmb->ncells1);
  u1.NewAthenaArray(NHYDRO,pmb->ncells3,pmb->ncells2,pmb->ncells1);
  u0.ZeroClear(); u1.ZeroClear();
  CopyActive(pmb,pmb->phydro->u,u0);
  CopyActive(pmb,pmb->phydro->u,u1);

  pmb->pcoord->AddCoordTermsDivergence(dt,pmb->phydro->flux,
                                       pmb->phydro->w,bb_cc,u1);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      const Real th=pmb->pcoord->x2v(j);
      const Real cot=std::cos(th)/std::sin(th);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real r=pmb->pcoord->x1v(i);
        Check("flat-static","dD",k,j,i,
              u1(IDN,k,j,i)-u0(IDN,k,j,i),0.0,tol,max_err,failed,log);
        Check("flat-static","dTau",k,j,i,
              u1(IEN,k,j,i)-u0(IEN,k,j,i),0.0,tol,max_err,failed,log);
        Check("flat-static","dS1",k,j,i,
              u1(IM1,k,j,i)-u0(IM1,k,j,i),dt*2.0*p/r,
              tol,max_err,failed,log);
        Check("flat-static","dS2",k,j,i,
              u1(IM2,k,j,i)-u0(IM2,k,j,i),dt*p*cot,
              tol,max_err,failed,log);
        Check("flat-static","dS3",k,j,i,
              u1(IM3,k,j,i)-u0(IM3,k,j,i),0.0,tol,max_err,failed,log);
      }
    }
  }
  std::cout << "  case: flat-static done" << std::endl;
}

void RunSchwarzschildMoving(MeshBlock *pmb, AthenaArray<Real> &bb_cc,
                            Real dt, Real tol, Real &max_err,
                            bool &failed, std::stringstream &log) {
  const Real gamma_adi=pmb->peos->GetGamma();
  const Real M=0.20;
  const Real rho=1.10;
  const Real p=0.13;
  const Real uu1=0.22;
  const Real uu2=-0.035;
  const Real uu3=0.028;

  SetMetricState(pmb,M,0.0,0.0);
  SetPrimitive(pmb,rho,p,uu1,uu2,uu3);
  P2C(pmb,bb_cc);

  AthenaArray<Real> u0,u1;
  u0.NewAthenaArray(NHYDRO,pmb->ncells3,pmb->ncells2,pmb->ncells1);
  u1.NewAthenaArray(NHYDRO,pmb->ncells3,pmb->ncells2,pmb->ncells1);
  u0.ZeroClear(); u1.ZeroClear();
  CopyActive(pmb,pmb->phydro->u,u0);
  CopyActive(pmb,pmb->phydro->u,u1);

  pmb->pcoord->AddCoordTermsDivergence(dt,pmb->phydro->flux,
                                       pmb->phydro->w,bb_cc,u1);

  const Real wgas=rho+gamma_adi/(gamma_adi-1.0)*p;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      const Real th=pmb->pcoord->x2v(j);
      const Real st=std::sin(th);
      const Real ct=std::cos(th);
      const Real st2=st*st;
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real r=pmb->pcoord->x1v(i);
        const Real r2=r*r;
        const Real f=1.0-2.0*M/r;
        const Real alpha=std::sqrt(f);
        const Real qs=1.0/alpha;
        const Real g11=1.0/f;
        const Real g22=r2;
        const Real g33=r2*st2;
        const Real gi00=-1.0/f;
        const Real gi11=f;
        const Real gi22=1.0/r2;
        const Real gi33=1.0/(r2*st2);
        const Real W=std::sqrt(1.0+g11*uu1*uu1+g22*uu2*uu2+g33*uu3*uu3);
        const Real u0c=W/alpha;

        const Real T00=wgas*u0c*u0c+p*gi00;
        const Real T11=wgas*uu1*uu1+p*gi11;
        const Real T22=wgas*uu2*uu2+p*gi22;
        const Real T33=wgas*uu3*uu3+p*gi33;
        const Real S1=alpha*wgas*u0c*uu1;

        const Real dg00=-2.0*M/r2;
        const Real dg11=-2.0*M/(r2*f*f);
        const Real dg22=2.0*r;
        const Real dg33=2.0*r*st2;
        const Real dth_g33=2.0*r2*st*ct;
        const Real dalpha=M/(alpha*r2);

        const Real src1=alpha*qs*0.5*(dg00*T00+dg11*T11+dg22*T22+dg33*T33);
        const Real src2=alpha*qs*0.5*dth_g33*T33;
        const Real srce=-qs*S1*dalpha;  // K_ij=0 in MonopoleGravityDriver

        Check("schwarzschild-moving","dD",k,j,i,
              u1(IDN,k,j,i)-u0(IDN,k,j,i),0.0,tol,max_err,failed,log);
        Check("schwarzschild-moving","dTau",k,j,i,
              u1(IEN,k,j,i)-u0(IEN,k,j,i),dt*srce,
              tol,max_err,failed,log);
        Check("schwarzschild-moving","dS1",k,j,i,
              u1(IM1,k,j,i)-u0(IM1,k,j,i),dt*src1,
              tol,max_err,failed,log);
        Check("schwarzschild-moving","dS2",k,j,i,
              u1(IM2,k,j,i)-u0(IM2,k,j,i),dt*src2,
              tol,max_err,failed,log);
        Check("schwarzschild-moving","dS3",k,j,i,
              u1(IM3,k,j,i)-u0(IM3,k,j,i),0.0,
              tol,max_err,failed,log);
      }
    }
  }
  std::cout << "  case: schwarzschild-moving done" << std::endl;
}

} // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  AllocateRealUserMeshDataField(2);
  ruser_mesh_data[0].NewAthenaArray(1);
  ruser_mesh_data[1].NewAthenaArray(1);
  ruser_mesh_data[0](0)=0.0;
  ruser_mesh_data[1](0)=0.0;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(ncells1+1);
  ruser_meshblock_data[1].NewAthenaArray(ncells1+1);
  ruser_meshblock_data[0].ZeroClear();
  ruser_meshblock_data[1].ZeroClear();
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real tol=pin->GetOrAddReal("problem","tol",2.0e-12);
  const Real dt_test=pin->GetOrAddReal("problem","dt_test",0.37);

  AthenaArray<Real> bb_cc;
  bb_cc.NewAthenaArray(3,ncells3,ncells2,ncells1);
  bb_cc.ZeroClear();

  if (gid==0) {
    Real max_err=0.0;
    bool failed=false;
    std::stringstream log;

    RunFlatStatic(this,bb_cc,dt_test,tol,max_err,failed,log);
    RunSchwarzschildMoving(this,bb_cc,dt_test,tol,max_err,failed,log);

    std::cout << std::setprecision(16)
              << "Valencia source unit test: max scaled error = "
              << max_err << std::endl;
    if (failed) {
      std::stringstream msg;
      msg << "### FATAL ERROR: Valencia source unit test failed\n" << log.str();
      ATHENA_ERROR(msg);
    }
    std::cout << "Valencia source unit test: PASS" << std::endl;
  }

  SetMetricState(this,0.0,0.0,0.0);
  SetPrimitive(this,1.0,0.1,0.0,0.0,0.0);
  P2C(this,bb_cc);
}

void MeshBlock::InitializeAtmosphere(ParameterInput *pin) {
  return;
}

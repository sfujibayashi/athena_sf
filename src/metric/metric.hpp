#ifndef METRIC_METRIC_HPP_
#define METRIC_METRIC_HPP_

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"

class MeshBlock;
class ParameterInput;

class Metric {
private:
  Real bh_mass_, bh_spin_;
  void InvertSpatialMetric(Real g11, Real g12, Real g13,
                           Real g22, Real g23, Real g33,
                           Real &gi11, Real &gi12, Real &gi13,
                           Real &gi22, Real &gi23, Real &gi33) const ;

  Real DetSpatialMetric(Real g11, Real g12, Real g13,
                        Real g22, Real g23, Real g33) const ;
  
  void ConstructCovariantMetric(
    Real r, Real theta, Real phi, Real Psi, Real dm,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const;

  void ConstructCellCovariantMetric(
    int k, int j, int i,
    Real &g00, Real &g01, Real &g02, Real &g03,
    Real &g11, Real &g12, Real &g13,
    Real &g22, Real &g23, Real &g33) const;

  void InvertMetric(
    int i,
    const AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv) const;

public:

  Metric(MeshBlock *pmb, ParameterInput *pin);
  ~Metric();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Field

  // Quantities needed to construct perturbed field (Only monopole l=0 mode)
  AthenaArray<Real> Psi_;
  AthenaArray<Real> delta_m_;

  AthenaArray<Real> Psi_face1_;
  AthenaArray<Real> delta_m_face1_;

  void Update(Real time);

  void CellMetric(const int k, const int j,
                  const int il, const int iu,
                  AthenaArray<Real> &g,
                  AthenaArray<Real> &g_inv);

  void Face1Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  
  Real SqrtMinusG(int k, int j, int i) const;

  Real DensitizationFactor(int k, int j, int i) const;

  void SetBlackHoleMass(Real mass);

  Real GetBlackHoleMass() const;


};

#endif  // METRIC_METRIC_HPP_

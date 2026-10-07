#ifndef METRIC_METRIC_HPP_
#define METRIC_METRIC_HPP_

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"

class MeshBlock;
class ParameterInput;

class Metric {
private:
  
  void InvertSpatialMetric(Real g11, Real g12, Real g13,
                           Real g22, Real g23, Real g33,
                           Real &gi11, Real &gi12, Real &gi13,
                           Real &gi22, Real &gi23, Real &gi33) const ;

  Real DetSpatialMetric(Real g11, Real g12, Real g13,
                        Real g22, Real g23, Real g33) const ;
  
public:

  Metric(MeshBlock *pmb, ParameterInput *pin);
  ~Metric();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Field

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

  Real CellDensitizationFactor(int k, int j, int i) const;
  Real Face1DensitizationFactor(int k, int j, int i) const;
  Real Face2DensitizationFactor(int k, int j, int i) const;
  Real Face3DensitizationFactor(int k, int j, int i) const;

  Real CellSpatialDensitizationFactor(int k, int j, int i) const;
  Real Face1SpatialDensitizationFactor(int k, int j, int i) const;
  Real Face2SpatialDensitizationFactor(int k, int j, int i) const;
  Real Face3SpatialDensitizationFactor(int k, int j, int i) const;

  void CellExtrinsicCurvature(const int k, const int j,
                  const int il, const int iu,
                  AthenaArray<Real> &k11, AthenaArray<Real> &k22, AthenaArray<Real> &k33) const;
  
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


  void CellMetricRadialDerivatives(const int k, const int j,
				   const int il, const int iu,
				   AthenaArray<Real> &d1_g00,
				   AthenaArray<Real> &d1_g11) const;


};

#endif  // METRIC_METRIC_HPP_

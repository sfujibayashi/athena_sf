#ifndef METRIC_METRIC_HPP_
#define METRIC_METRIC_HPP_

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"

class MeshBlock;
class ParameterInput;

class Metric {
private:
  Real bh_mass_;
  void InvertSpatialMetric(Real g11, Real g12, Real g13,
                           Real g22, Real g23, Real g33,
                           Real &gi11, Real &gi12, Real &gi13,
                           Real &gi22, Real &gi23, Real &gi33);
  
public:
  
  enum {
    I_G11 = 0,
    I_G12,
    I_G13,
    I_G22,
    I_G23,
    I_G33,
    N_GAMMA
  };

  Metric(MeshBlock *pmb, ParameterInput *pin);
  ~Metric();

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Field

  // ADM metric variables
  AthenaArray<Real> alpha;
  AthenaArray<Real> beta;
  AthenaArray<Real> gamma;

  void Update(Real time);

  void CellMetric(const int k, const int j,
                  const int il, const int iu,
                  AthenaArray<Real> &g,
                  AthenaArray<Real> &g_inv);
  
};

#endif  // METRIC_METRIC_HPP_

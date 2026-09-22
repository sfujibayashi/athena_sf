#ifndef METRIC_METRIC_HPP_
#define METRIC_METRIC_HPP_

#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParameterInput;

class Metric {
 public:

  enum {
    I_GXX = 0,
    I_GXY,
    I_GXZ,
    I_GYY,
    I_GYZ,
    I_GZZ,
    N_GAMMA
  };

  Metric(MeshBlock *pmb, ParameterInput *pin);
  ~Metric();

  void Update(Real time);

  MeshBlock *pmy_block;  // ptr to MeshBlock containing this Field

  // ADM metric variables
  AthenaArray<Real> alpha;
  AthenaArray<Real> beta;
  AthenaArray<Real> gamma;
};

#endif  // METRIC_METRIC_HPP_

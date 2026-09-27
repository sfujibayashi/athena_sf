#include "monopole_gravity.hpp"

#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../metric/metric.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"

MonopoleGravity::MonopoleGravity(Mesh *pm, ParameterInput *pin)
  : pmy_mesh_(pm) {

  nr_ = pm->mesh_size.nx1;

  dm_shell_global_.NewAthenaArray(nr_);
  delta_m_face_global_.NewAthenaArray(nr_ + 1);
  psi_face_global_.NewAthenaArray(nr_ + 1);
  
  dm_shell_global_.ZeroClear();
  delta_m_face_global_.ZeroClear();
  psi_face_global_.ZeroClear();
}

MonopoleGravity::~MonopoleGravity() {
  dm_shell_global_.DeleteAthenaArray();
  delta_m_face_global_.DeleteAthenaArray();
  psi_face_global_.DeleteAthenaArray();
}

void MonopoleGravity::Update(){
}

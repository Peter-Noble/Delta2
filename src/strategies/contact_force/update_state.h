#include "../../model/state.h"
#include "../../model/particle.h"

namespace Delta2 {
    Delta2::State updateState(Delta2::State& future, double t, const Eigen::Vector3d& impulse, const Eigen::Vector3d& impulse_offset, const Eigen::Vector3d& rotational_impulse, const Eigen::Vector3d& rotational_impulse_offset, Delta2::Particle* p);
}

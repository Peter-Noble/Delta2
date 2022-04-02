#pragma once

#include "particle.h"

namespace Delta2 {
    namespace model {
        void staticApply(Delta2::Particle& P, double t, bool apply=true);
        double integrateEuler(Delta2::Particle& P, const Eigen::Vector3d& force, const Eigen::Vector3d& torque, double t, bool apply=true);
    }
}

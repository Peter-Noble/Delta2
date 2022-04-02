#pragma once
#include "particle.h"

namespace Delta2 {
    namespace model {
        Eigen::Matrix3d tetrahedronInertiaTensor(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d);

        Eigen::Matrix3d pointInertiaTensor(const Eigen::Vector3d& pt);
    }
}

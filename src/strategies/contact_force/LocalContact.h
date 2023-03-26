#pragma once

#include <Eigen/Dense>

namespace Delta2 {
    struct LocalContact {
        Eigen::Vector3d local_A;
        Eigen::Vector3d local_B;
        Eigen::Vector3d global_centre;
        Eigen::Vector3d global_normal;
        Eigen::Vector3d global_tangent;
        double mass_eff_normal;
        double mass_eff_tangent;

        double start_dist;

        double impulse;
        double impulse_offset;

        Eigen::Vector3d last_friction_impulse;
        Eigen::Vector3d last_friction_rotational_impulse_A;
        Eigen::Vector3d last_friction_rotational_impulse_B;
        Eigen::Vector3d last_friction_impulse_fd;
        Eigen::Vector3d last_friction_rotational_impulse_A_fd;
        Eigen::Vector3d last_friction_rotational_impulse_B_fd;
    };
}

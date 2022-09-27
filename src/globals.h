#pragma once

#include <vector>
#include <Eigen/Dense>
#include <mutex>

namespace Delta2 {
    namespace globals {
        extern std::mutex contact_draws_lock;
        extern std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> contact_draws;
    }
}

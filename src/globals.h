#pragma once

#include <vector>
#include <Eigen/Dense>
#include <mutex>
#include "./common/cli_options.h"

namespace Delta2 {
    namespace globals {
        extern std::mutex contact_draws_lock;
        extern std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> contact_draws;
        extern Delta2::common::Options opt;
    }
}

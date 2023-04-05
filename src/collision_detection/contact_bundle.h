#pragma once

#include <Eigen/Dense>
#include <vector>
#include "contact.h"

namespace Delta2 {
    namespace collision {
        struct ContactBundle {
            uint32_t lower;
            uint32_t upper;
            std::vector<int> hits;
            Eigen::Vector3d normal_average;
            Eigen::Vector3d tangent_average;
        };
    }
}

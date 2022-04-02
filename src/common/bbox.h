#pragma once

#include <Eigen/Dense>

namespace Delta2 {
    namespace common {
        template<typename real>
        bool overlappingRange(real a_lower, real a_upper, real b_lower, real b_upper) {
            if (b_lower < a_lower) {
                if (b_upper > a_lower) {
                    return true;
                }
                else {
                    return false;
                }
            }
            else if (b_lower < a_upper) {
                return true;
            }
            else {
                return false;
            }
        }

        template<typename real>
        struct bbox {
            Eigen::Vector<real, 3> lower;
            Eigen::Vector<real, 3> upper;
            bool overlap(bbox<real> other) {
                bool x = overlappingRange(lower.x(), upper.x(), other.lower.x(), other.upper.x());
                bool y = overlappingRange(lower.y(), upper.y(), other.lower.y(), other.upper.y());
                bool z = overlappingRange(lower.z(), upper.z(), other.lower.z(), other.upper.z());
                return x && y && z;
            }
            bbox<real> expand(bbox<real> other) {
                bbox<real> result;
                result.lower = {std::min(lower.x(), other.lower.x()), std::min(lower.y(), other.lower.y()), std::min(lower.z(), other.lower.z())};
                result.upper = {std::max(upper.x(), other.upper.x()), std::max(upper.y(), other.upper.y()), std::max(upper.z(), other.upper.z())};
                return result;
            }
            bbox<real> expand(real eps) {
                bbox<real> result;
                result.lower = {lower.x() - eps, lower.y() - eps, lower.z() - eps};
                result.upper = {upper.x() + eps, upper.y() + eps, upper.z() + eps};
                return result;
            }
        };
    }
}

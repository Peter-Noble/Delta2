#include <vector>
#include "triangle.h"

namespace Delta2 {
    namespace common {
        std::vector<Triangle<double>> convex_hull(std::vector<Eigen::Vector3d> pts);
    }
}

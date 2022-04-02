#include "convex_hull.h"
#include "../quickhull/QuickHull.hpp"

namespace Delta2 {
    namespace common {
        std::vector<Triangle<double>> convex_hull(std::vector<Eigen::Vector3d> pts) {
            // using namespace quickhull;
            quickhull::QuickHull<double> qh; // Could be double as well
            std::vector<quickhull::Vector3<double>> pointCloud;
            // Add points to point cloud
            for (const Eigen::Vector3d& p : pts) {
                pointCloud.push_back(quickhull::Vector3<double>(p.x(), p.y(), p.z()));
            }

            auto hull = qh.getConvexHull(pointCloud, true, true);
            const std::vector<size_t>& indexBuffer = hull.getIndexBuffer();
            //const auto& vertexBuffer = hull.getVertexBuffer();
            // Do what you want with the convex triangle mesh

            std::vector<Triangle<double>> result;
            for (int t_i = 0; t_i < indexBuffer.size(); t_i+=3) {
                Triangle<double> t(pts[t_i], pts[t_i+1], pts[t_i+2]);
                result.push_back(t);
            }

            return result;
        }
    }
}

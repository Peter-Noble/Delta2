#pragma once

#include <Eigen/Dense>

#include "triangle.h"

namespace Delta2 {
    namespace common {
        template<typename real, int vertices, int faces, int dimensions, int ngon>
        std::vector<Triangle<real>> toTriangles(const Eigen::Matrix<real, vertices, dimensions>& V, const Eigen::Matrix<int, faces, ngon>& F, int tris) {
            std::vector<Triangle<real>> result;
            for (int r = 0; r < tris; r++) {
                Eigen::Vector<real, 3> A(V(F(r,0),0), V(F(r,0),1), V(F(r,0),2));
                Eigen::Vector<real, 3> B(V(F(r,1),0), V(F(r,1),1), V(F(r,1),2));
                Eigen::Vector<real, 3> C(V(F(r,2),0), V(F(r,2),1), V(F(r,2),2));
                Triangle<real> T(A, B, C);
                result.push_back(T);
            }
            return result;
        }

        template<typename real>
        bool isCoplanar(const Eigen::Matrix<real, -1, -1> V, double tol = 1e-8) {
            if (V.size() < 4) {
                return true;
            }
            Delta2::common::Triangle<real> T({V(0,0), V(0,1), V(0,2)}, {V(1,0), V(1,1), V(1,2)}, {V(2,0), V(2,1), V(2,2)});

            for (int r = 3; r < V.rows(); r++) {
                Eigen::Vector<real, 3> pt({V(r,0), V(r,1), V(r,2)});
                bool _;
                Eigen::Vector<real, 3> proj = T.projectToPlane(pt, _);
                if ((pt - proj).norm() > tol) {
                    return false;
                }
            }
            return true;
        }
    }
}

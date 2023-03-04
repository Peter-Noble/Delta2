#if !defined(NOGL)
#pragma once

#include "igl/opengl/glfw/Viewer.h"
#include "igl/PI.h"
#include "../common/triangle.h"
#include "../common/basic_utils.h"

namespace Delta2 {
    namespace debug {
        template<typename real>
        void debugViewPoints(std::vector<Eigen::Vector<real, 3>> pts) {
            Eigen::Matrix<real, -1, -1> V;
            V.resize(pts.size(), 3);
            for (int i = 0; i < pts.size(); i++) {
                V.row(i) << pts[i].x(), pts[i].y(), pts[i].z();
            }

            Eigen::Matrix<real, 4, 4> T = Delta2::common::transformationMatrix(Eigen::Vector<real, 3>({-igl::PI/2.0, 0.0, 0.0}), Eigen::Vector<real, 3>({0.0, 0.0, 0.0}));
            V = common::transform(V, T);

            Eigen::MatrixXd C;
            C.resize(1, 3);
            C << 1.0, 0.0, 0.0;

            igl::opengl::glfw::Viewer viewer;
            viewer.data().add_points(V.template cast<double>(), C);

            viewer.launch(true, false, "Delta2 - debug", 0, 0);
        }
    }
}
#endif
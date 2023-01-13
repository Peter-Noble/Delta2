#pragma once

#include <Eigen/Dense>

namespace Delta2 {
    namespace common {
        void cube(double r, Eigen::Matrix<double, -1, -1>& V, Eigen::MatrixXi& F);
        void sphere(const double& r, const int res, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
        void plane(int width, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    }
}
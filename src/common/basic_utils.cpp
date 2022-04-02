#include "basic_utils.h"

void Delta2::common::concatMesh(const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA, const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB, Eigen::MatrixXd& V_out, Eigen::MatrixXi& F_out) {
    if (VA.rows() == 0) {
        V_out = VB;
        F_out = FB;
    }
    else if (VB.rows() == 0) {
        V_out = VA;
        F_out = FB;
    }
    else {
        Eigen::MatrixXd V(VA.rows() + VB.rows(), VA.cols());
        V << VA, VB;
        Eigen::MatrixXi F(FA.rows() + FB.rows(), FA.cols());
        F << FA, (FB.array() + VA.rows());
        V_out = V;
        F_out = F;
    }
}

void Delta2::common::concatColors(const Eigen::MatrixXd& CA, const Eigen::MatrixXd& CB, Eigen::MatrixXd& C_out) {
    if (CA.rows() == 0) {
        C_out = CB;
    }
    else if (CB.rows() == 0) {
        C_out = CA;
    }
    else {
        Eigen::MatrixXd C(CA.rows() + CB.rows(), CA.cols());
        C << CA, CB;
        C_out = C;
    }
}

Eigen::Quaterniond Delta2::common::exp(Eigen::Vector3d d) {
    // https://math.stackexchange.com/questions/400028/quaternion-exponential-map-rotations-and-interpolation
    if (d.norm() == 0.0) {
        return Eigen::Quaterniond::Identity();
    }
    Eigen::Vector3d v = d / d.norm() * std::sin(d.norm());
    Eigen::Quaterniond result = Eigen::Quaterniond(std::cos(d.norm()), v.x(), v.y(), v.z());
    return result;
}

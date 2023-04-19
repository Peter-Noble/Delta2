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

// float fast_sin(float x)
// {
//     const double PI = 3.14159265358979323846264338327950288;
//     const float B = 4 / PI;
//     const float C = -4 / (PI*PI);

//     return -(B * x + C * x * ((x < 0) ? -x : x));
// } 

// Eigen::Quaterniond Delta2::common::exp(Eigen::Vector3d d) {
//     // https://math.stackexchange.com/questions/400028/quaternion-exponential-map-rotations-and-interpolation
//     double norm = d.norm();
//     if (norm == 0.0) {
//         return Eigen::Quaterniond::Identity();
//     }
//     // double sin = norm < 0.01 ? norm : std::sin(d.norm());  // Small angle approximation.  Doesn't really make a big difference
//     double sin = std::sin(d.norm());
//     // double sin = fast_sin(d.norm());
//     Eigen::Vector3d v = d / norm * sin;
//     Eigen::Quaterniond result = Eigen::Quaterniond(1.0-sin, v.x(), v.y(), v.z());
//     return result;
// }

Eigen::Quaterniond Delta2::common::exp(Eigen::Vector3d d) {
    // https://math.stackexchange.com/questions/400028/quaternion-exponential-map-rotations-and-interpolation
    if (d.norm() == 0.0) {
        return Eigen::Quaterniond::Identity();
    }
    Eigen::Vector3d v = d / d.norm() * std::sin(d.norm());
    Eigen::Quaterniond result = Eigen::Quaterniond(std::cos(d.norm()), v.x(), v.y(), v.z());
    return result;
}
// Things that only depend on basic data types
// ie. not allowed to include other header files from Delta
#pragma once
#include <Eigen/Dense>
#include <iostream>
#include "bbox.h"

namespace Delta2 {
    namespace common {
        template<typename real>
        real clamp01(real t) {
            if (t < 0) {
                return 0;
            } else if (t > 1) {
                return 1;
            } else {
                return t;
            }
        }
        template float clamp01(float);
        template double clamp01(double);

        template<typename real>
        Eigen::Vector<real, 3> lerp(Eigen::Vector<real, 3> A, Eigen::Vector<real, 3> B, real t) {
            return A * (1 - t) + B * t;
        }
        template Eigen::Vector3f lerp(Eigen::Vector3f, Eigen::Vector3f, float);
        template Eigen::Vector3d lerp(Eigen::Vector3d, Eigen::Vector3d, double);

        template<typename real>
        real lerp(real A, real B, real t) {
            return A * (1 - t) + B * t;
        }
        template float lerp(float, float, float);
        template double lerp(double, double, double);

        template<typename real>
        Eigen::Vector<real, 3> transform(const Eigen::Vector<real, 3>& pt, const Eigen::Matrix<real, 4, 4>& trans) {
            Eigen::Vector<real, 4> tmp = trans * pt.homogeneous();
            Eigen::Vector<real, 3> result = tmp.template head<3>() / tmp.w();
            return result;
        }
        template Eigen::Vector3f transform(const Eigen::Vector3f&, const Eigen::Matrix4f&);
        template Eigen::Vector3d transform(const Eigen::Vector3d&, const Eigen::Matrix4d&);

        template<typename real>
        Eigen::Vector<real, 3> transform(const Eigen::Vector<real, 3>& pt, const Eigen::Matrix<real, 3, 3>& trans) {
            Eigen::Vector<real, 3> result = trans * pt;
            return result;
        }
        template Eigen::Vector3f transform(const Eigen::Vector3f&, const Eigen::Matrix3f&);
        template Eigen::Vector3d transform(const Eigen::Vector3d&, const Eigen::Matrix3d&);

        template<typename real>
        Eigen::Matrix<real, -1, -1> transform(const Eigen::Matrix<real, -1, -1>& V, const Eigen::Matrix<real, 4, 4>& trans) {
            Eigen::Matrix<real, -1, -1> result;
            result.resize(V.rows(), 3);
            for (int i = 0; i < V.rows(); i++) {
                Eigen::Vector<real, 3> t = transform(Eigen::Vector<real, 3>({V(i, 0), V(i, 1), V(i, 2)}), trans);
                result(i, 0) = t.x();
                result(i, 1) = t.y();
                result(i, 2) = t.z();
            }
            return result;
        }
        template Eigen::Matrix<float, -1, -1> transform(const Eigen::Matrix<float, -1, -1>& V, const Eigen::Matrix<float, 4, 4>& trans);
        template Eigen::Matrix<double, -1, -1> transform(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<double, 4, 4>& trans);


        template<typename real>
        Eigen::Quaternion<real> eulerAnglesToQuaternion(const Eigen::Vector<real, 3> angles) {
            Eigen::Quaternion<real> qa;
            qa = Eigen::AngleAxis<real>(angles[0], Eigen::Vector<real, 3>::UnitX())
                * Eigen::AngleAxis<real>(angles[1], Eigen::Vector<real, 3>::UnitY())
                * Eigen::AngleAxis<real>(angles[2], Eigen::Vector<real, 3>::UnitZ());
            return qa;
        }

        template<typename real>
        Eigen::Matrix<real, 4, 4> transformationMatrix(const Eigen::Quaternion<real>& rotation, const Eigen::Vector<real, 3>& translation) {
            Eigen::Matrix<real, 4, 4> transform;
            transform.setIdentity();
            transform.template block<3,3>(0,0) = rotation.toRotationMatrix();
            transform.template block<3,1>(0,3) = translation;
            return transform;
        }
        template Eigen::Matrix4f transformationMatrix(const Eigen::Quaternionf&, const Eigen::Vector3f&);
        template Eigen::Matrix4d transformationMatrix(const Eigen::Quaterniond&, const Eigen::Vector3d&);

        template<typename real>
        Eigen::Matrix<real, 4, 4> transformationMatrix(const Eigen::Vector<real, 3>& rotation, const Eigen::Vector<real, 3>& translation) {
            return transformationMatrix(eulerAnglesToQuaternion(rotation), translation);
        }
        template Eigen::Matrix4f transformationMatrix(const Eigen::Vector3f&, const Eigen::Vector3f&);
        template Eigen::Matrix4d transformationMatrix(const Eigen::Vector3d&, const Eigen::Vector3d&);

        void concatMesh(const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA, const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB, Eigen::MatrixXd& V_out, Eigen::MatrixXi& F_out);
        void concatColors(const Eigen::MatrixXd& CA, const Eigen::MatrixXd& CB, Eigen::MatrixXd& C_out);
        Eigen::Quaterniond exp(Eigen::Vector3d d);

        template<typename real>
        real degToRad(real deg) {
            real pi = 3.14159265358979323846;
            return 2.0 * pi * deg / 360.0;
        }
    }
}
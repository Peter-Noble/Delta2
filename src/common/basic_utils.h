// Things that only depend on basic data types
// ie. not allowed to include other header files from Delta
#pragma once
#include <Eigen/Dense>
#include <iostream>
#include "bbox.h"

namespace Delta2 {
    namespace common {
        inline float clamp01(float t) {
            if (t < 0) {
                return 0;
            } else if (t > 1) {
                return 1;
            } else {
                return t;
            }
        }
        inline double clamp01(double t) {
            if (t < 0) {
                return 0;
            } else if (t > 1) {
                return 1;
            } else {
                return t;
            }
        }

        inline Eigen::Vector<float, 3> lerp(const Eigen::Vector<float, 3>& A, const Eigen::Vector<float, 3>& B, float t) {
            return A * (1 - t) + B * t;
        }
        inline Eigen::Vector<double, 3> lerp(const Eigen::Vector<double, 3>& A, const Eigen::Vector<double, 3>& B, double t) {
            return A * (1 - t) + B * t;
        }
        
        inline float lerp(float A, float B, float t) {
            return A * (1 - t) + B * t;
        }
        inline double lerp(double A, double B, double t) {
            return A * (1 - t) + B * t;
        }

        inline Eigen::Vector<float, 3> transform(const Eigen::Vector<float, 3>& pt, const Eigen::Matrix<float, 4, 4>& trans) {
            Eigen::Vector<float, 4> tmp = trans * pt.homogeneous();
            Eigen::Vector<float, 3> result = tmp.template head<3>() / tmp.w();
            return result;
        }
        inline Eigen::Vector<double, 3> transform(const Eigen::Vector<double, 3>& pt, const Eigen::Matrix<double, 4, 4>& trans) {
            Eigen::Vector<double, 4> tmp = trans * pt.homogeneous();
            Eigen::Vector<double, 3> result = tmp.template head<3>() / tmp.w();
            return result;
        }

        inline Eigen::Vector<float, 3> transform(const Eigen::Vector<float, 3>& pt, const Eigen::Matrix<float, 3, 3>& trans) {
            Eigen::Vector<float, 3> result = trans * pt;
            return result;
        }
        inline Eigen::Vector<double, 3> transform(const Eigen::Vector<double, 3>& pt, const Eigen::Matrix<double, 3, 3>& trans) {
            Eigen::Vector<double, 3> result = trans * pt;
            return result;
        }

        inline Eigen::Vector<float, 3> transform(const Eigen::Vector<float, 3>& pt, const Eigen::Quaternion<float>& trans) {
            Eigen::Vector<float, 3> result = trans.toRotationMatrix() * pt;
            return result;
        }
        inline Eigen::Vector<double, 3> transform(const Eigen::Vector<double, 3>& pt, const Eigen::Quaternion<double>& trans) {
            Eigen::Vector<double, 3> result = trans.toRotationMatrix() * pt;
            return result;
        }

        inline Eigen::Matrix<float, -1, -1> transform(const Eigen::Matrix<float, -1, -1>& V, const Eigen::Matrix<float, 4, 4>& trans) {
            Eigen::Matrix<float, -1, -1> result;
            result.resize(V.rows(), 3);
            for (int i = 0; i < V.rows(); i++) {
                Eigen::Vector<float, 3> t = transform(Eigen::Vector<float, 3>({V(i, 0), V(i, 1), V(i, 2)}), trans);
                result(i, 0) = t.x();
                result(i, 1) = t.y();
                result(i, 2) = t.z();
            }
            return result;
        }
        inline Eigen::Matrix<double, -1, -1> transform(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<double, 4, 4>& trans) {
            Eigen::Matrix<double, -1, -1> result;
            result.resize(V.rows(), 3);
            for (int i = 0; i < V.rows(); i++) {
                Eigen::Vector<double, 3> t = transform(Eigen::Vector<double, 3>({V(i, 0), V(i, 1), V(i, 2)}), trans);
                result(i, 0) = t.x();
                result(i, 1) = t.y();
                result(i, 2) = t.z();
            }
            return result;
        }

        inline Eigen::Matrix<float, -1, -1> transform(const Eigen::Matrix<float, -1, -1>& V, const Eigen::Quaternion<float>& trans) {
            Eigen::Matrix<float, -1, -1> result;
            result.resize(V.rows(), 3);
            for (int i = 0; i < V.rows(); i++) {
                Eigen::Vector<float, 3> t = transform(Eigen::Vector<float, 3>({V(i, 0), V(i, 1), V(i, 2)}), trans.toRotationMatrix());
                result(i, 0) = t.x();
                result(i, 1) = t.y();
                result(i, 2) = t.z();
            }
            return result;
        }
        inline Eigen::Matrix<double, -1, -1> transform(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Quaternion<double>& trans) {
            Eigen::Matrix<double, -1, -1> result;
            result.resize(V.rows(), 3);
            for (int i = 0; i < V.rows(); i++) {
                Eigen::Vector<double, 3> t = transform(Eigen::Vector<double, 3>({V(i, 0), V(i, 1), V(i, 2)}), trans.toRotationMatrix());
                result(i, 0) = t.x();
                result(i, 1) = t.y();
                result(i, 2) = t.z();
            }
            return result;
        }

        inline Eigen::Quaternion<float> eulerAnglesToQuaternion(const Eigen::Vector<float, 3> angles) {
            Eigen::Quaternion<float> qa;
            qa = Eigen::AngleAxis<float>(angles[0], Eigen::Vector<float, 3>::UnitX())
                * Eigen::AngleAxis<float>(angles[1], Eigen::Vector<float, 3>::UnitY())
                * Eigen::AngleAxis<float>(angles[2], Eigen::Vector<float, 3>::UnitZ());
            return qa;
        }
        inline Eigen::Quaternion<double> eulerAnglesToQuaternion(const Eigen::Vector<double, 3> angles) {
            Eigen::Quaternion<double> qa;
            qa = Eigen::AngleAxis<double>(angles[0], Eigen::Vector<double, 3>::UnitX())
                * Eigen::AngleAxis<double>(angles[1], Eigen::Vector<double, 3>::UnitY())
                * Eigen::AngleAxis<double>(angles[2], Eigen::Vector<double, 3>::UnitZ());
            return qa;
        }

        inline Eigen::Matrix<float, 4, 4> transformationMatrix(const Eigen::Quaternion<float>& rotation, const Eigen::Vector<float, 3>& translation) {
            Eigen::Matrix<float, 4, 4> transform;
            transform.setIdentity();
            transform.template block<3,3>(0,0) = rotation.toRotationMatrix();
            transform.template block<3,1>(0,3) = translation;
            return transform;
        }
        inline Eigen::Matrix<double, 4, 4> transformationMatrix(const Eigen::Quaternion<double>& rotation, const Eigen::Vector<double, 3>& translation) {
            Eigen::Matrix<double, 4, 4> transform;
            transform.setIdentity();
            transform.template block<3,3>(0,0) = rotation.toRotationMatrix();
            transform.template block<3,1>(0,3) = translation;
            return transform;
        }
        
        inline Eigen::Matrix<float, 4, 4> transformationMatrix(const Eigen::Vector<float, 3>& rotation, const Eigen::Vector<float, 3>& translation) {
            return transformationMatrix(eulerAnglesToQuaternion(rotation), translation);
        }
        inline Eigen::Matrix<double, 4, 4> transformationMatrix(const Eigen::Vector<double, 3>& rotation, const Eigen::Vector<double, 3>& translation) {
            return transformationMatrix(eulerAnglesToQuaternion(rotation), translation);
        }
        
        void concatMesh(const Eigen::MatrixXd& VA, const Eigen::MatrixXi& FA, const Eigen::MatrixXd& VB, const Eigen::MatrixXi& FB, Eigen::MatrixXd& V_out, Eigen::MatrixXi& F_out);
        void concatColors(const Eigen::MatrixXd& CA, const Eigen::MatrixXd& CB, Eigen::MatrixXd& C_out);
        Eigen::Quaterniond exp(Eigen::Vector3d d);

        inline float degToRad(float deg) {
            float pi = 3.14159265358979323846;
            return 2.0 * pi * deg / 360.0;
        }
        inline double degToRad(double deg) {
            double pi = 3.14159265358979323846;
            return 2.0 * pi * deg / 360.0;
        }
    }
}
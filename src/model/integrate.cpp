#include "integrate.h"
#include "inertia.h"
#include "../common/basic_utils.h"

namespace Delta2 {
    namespace model {
        void staticApply(Delta2::Particle& P, double t, bool apply) {
            if (apply) {
                P.current_state.setTime(P.current_state.getTime() + t);
                P.future_state = P.current_state;
            } else {
                P.future_state = P.current_state;
                P.future_state.setTime(P.current_state.getTime() + t);
            }
        }
        double integrateEuler(Delta2::Particle& P, const Eigen::Vector3d& force, const Eigen::Vector3d& torque, double t, bool apply) {
            Delta2::State updated;
            updated.setTime(P.current_state.getTime() + t);
            updated.setVelocity(P.current_state.getVelocity() + t * force / P.getMass());
            updated.setTranslation(P.current_state.getTranslation() + t * updated.getVelocity());

            // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

            Eigen::Vector3d angular_momentum = P.current_state.getAngularMomentum() + t * torque;

            Eigen::Quaterniond rotation = P.current_state.getRotation();

            updated.setAngular(angular_momentum);

            Eigen::Matrix3d R = rotation.toRotationMatrix();
            Eigen::Matrix3d Iinv = R * P.getInverseInertiaMatrix() * R.transpose();
            Eigen::Vector3d omega = Iinv * angular_momentum;

            /*rotation.inverse() * P.getInverseInertiaMatrix() * rotation.toRotationMatrix() * angular_momentum

            Eigen::Quaterniond torque_global;
            torque_global.w() = 0.0;
            torque_global.x() = angular_momentum.x();
            torque_global.y() = angular_momentum.y();
            torque_global.z() = angular_momentum.z();

            Eigen::Quaterniond torque_local = rotation.conjugate() * torque_global * rotation;
            P.getInverseInertiaMatrix() */

            // Eigen::Quaterniond q_omega;
            // q_omega.w() = 0;
            // q_omega.x() = 0.5 * omega[0];
            // q_omega.y() = 0.5 * omega[1];
            // q_omega.z() = 0.5 * omega[2];

            //Eigen::Quaterniond a_ang = q_omega * rotation;
            //Eigen::Quaterniond a_ang(star(omega) * R);

            // rotation.w() += t * a_ang.w();
            // rotation.x() += t * a_ang.x();
            // rotation.y() += t * a_ang.y();
            // rotation.z() += t * a_ang.z();
            // rotation.normalize();
            
            Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

            /*Eigen::Quaterniond rot_delta;
            rot_delta.w() = std::cos(omega.norm() * t / 2.0);
            double fac = std::sin(omega.norm() * t / 2.0) / omega.norm();
            rot_delta.x() = omega.x() * fac;
            rot_delta.y() = omega.y() * fac;
            rot_delta.z() = omega.z() * fac;*/

            // rotation = rotation * rot_delta;
            rotation = rot_delta * rotation;

            // Eigen::Vector3d rea = rotation.toRotationMatrix().eulerAngles(1,2,3);
            updated.setRotation(rotation);

            if (!updated.isValid()) {
                printf("Invalid state\n");
            }

            double diff = P.maxDifferenceToCurrent(updated);
            if (apply) {
                P.last_state = P.current_state;
                P.current_state = updated;
                P.projectFutureState(t);
            } else {
                P.future_state = updated;
            }
            return diff;
        }
    }
}

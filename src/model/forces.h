#pragma once
#include <Eigen/Dense>
#include <functional>

#include "particle.h"
#include "particle_handler.h"
#include "../collision_detection/contact.h"

namespace Delta2 {
    namespace model {
        // Using force = (outer - inner) / (dist - inner)^2 - (outer - inner) / (outer - inner)^2
        // Integral of force = (inner - outer) / (dist - inner) - dist / (inner - outer) + C
        template<typename real>
        Eigen::Vector<real, 3> calcImpulse(const collision::Contact<real>& hit, real time_step) {
            // Impulse on A from hit.

            Eigen::Vector<real, 3> hit_normal = hit.A - hit.B;

            real inner = hit.eps_inner_a + hit.eps_inner_b;
            real outer = hit.eps_a + hit.eps_b;
            real dist_end = hit_normal.norm();

            Eigen::Vector<real, 3> a_vel = hit.p_a->futurePointVelocity(hit.A);
            Eigen::Vector<real, 3> b_vel = hit.p_b->futurePointVelocity(hit.B);

            Eigen::Vector<real, 3> rel_vel = a_vel - b_vel;
            real closing_velocity = std::max(0.0, -(rel_vel.dot(hit_normal.normalized())));

            real dist_start = std::min(dist_end - closing_velocity * time_step, outer);

            auto integral = [=](real dist) { return (inner - outer) / (dist - inner) - dist / (inner - outer); };

            real impulse = integral(dist_end) - integral(dist_start);

            real mass_weight;
            if (hit.p_a->is_static || hit.p_b->is_static) {
                mass_weight = 2.0;
            } else {
                real mass_ratio = std::max(hit.p_a->getMass(), hit.p_b->getMass()) / std::min(hit.p_a->getMass(), hit.p_b->getMass());
                mass_weight = (mass_ratio - 1.0) / mass_ratio + 1.0;
            }

            return mass_weight * impulse * hit_normal / hit_normal.norm();
        }

        template<typename real>
        Eigen::Vector<real, 3> calcForce(const collision::Contact<real>& hit, real multiplier=1.0, real time_step=1.0) {
            return calcImpulse(hit, time_step) / time_step * multiplier;
        }

        /*template<typename real>
        Eigen::Vector<real, 3> calcForce(const collision::Contact<real>& hit, real multiplier=1.0, real time_step=1.0) {
            // Force on A from hit.

            Eigen::Vector<real, 3> hit_normal = hit.A - hit.B;
            // printf("Hit dist: %f\n", hit_normal.norm());

            real inner = hit.eps_inner_a + hit.eps_inner_b;
            real outer = hit.eps_a + hit.eps_b;
            real dist = hit_normal.norm();

            real force = ((outer - inner) / (dist - inner) - 1.0) * multiplier;

            Eigen::Vector<real, 3> a_vel = hit.p_a->futurePointVelocity(hit.A);
            Eigen::Vector<real, 3> b_vel = hit.p_b->futurePointVelocity(hit.B);

            Eigen::Vector<real, 3> rel_vel = a_vel - b_vel;

            double closing_velocity = std::max(0.0, -(rel_vel.dot(hit_normal.normalized())));

            // double closing_velocity_proportion = closing_velocity * time_step / hit_normal.norm();
            double closing_velocity_proportion = closing_velocity / hit_normal.norm();

            double velocity_multiplier = closing_velocity_proportion + 1.0;

            double mass_multiplier = 0;
            if (!hit.p_a->is_static) {
                mass_multiplier += hit.p_a->getMass();
            }
            if (!hit.p_b->is_static) {
                mass_multiplier += hit.p_b->getMass();
            }
            if (hit.p_a->is_static || hit.p_b->is_static) {
                mass_multiplier *= 2;
            }

            if (rel_vel.dot(hit_normal) > 0.0) {
                force = 0.0;
                // printf("Rel vel cancel\n");
            }

            Eigen::Vector<real, 3> result = hit_normal.normalized() * force * velocity_multiplier * mass_multiplier;

            if (!result.allFinite()) {
                printf("Invalid force\n");
            }

            return result;
        }*/

        template<typename real>
        Eigen::Vector<real, 3> calcTorque(const Eigen::Vector<real, 3>& force, const Eigen::Vector<real, 3>& point, const Eigen::Vector<real, 3>& centre_of_mass) {
            //Eigen::Vector<real, 3> trans_force = -force - centre_of_mass;
            Eigen::Vector<real, 3> trans_force = -force;
            Eigen::Vector<real, 3> trans_point = point - centre_of_mass;
            Eigen::Vector<real, 3> torque = trans_force.cross(trans_point);
            return torque;
        }

        template<typename real>
        void calcFriction(Delta2::collision::Contact<real> hit, const Delta2::Particle& A, const Delta2::Particle& B, Eigen::Vector<real, 3>& a_force_out, Eigen::Vector<real, 3>& a_torque_out, Eigen::Vector<real, 3>& b_force_out, Eigen::Vector<real, 3>& b_torque_out) {
            if (std::abs(A.future_state.getTime() - B.future_state.getTime()) > 1e-8) {
                throw std::invalid_argument("Input particles must have futures states with the same time");
            }

            // Eigen::Vector<real, 3> a_rel_vel_new = A.current_state.pointVelocity(hit.A, A.future_state);
            // Eigen::Vector<real, 3> b_rel_vel_new = B.current_state.pointVelocity(hit.B, B.future_state);

            Eigen::Vector<real, 3> a_rel_vel_new = A.futurePointVelocity(hit.A);
            Eigen::Vector<real, 3> b_rel_vel_new = B.futurePointVelocity(hit.B);

            Eigen::Vector<real, 3> rel_vel_new = a_rel_vel_new - b_rel_vel_new;

            real mu = sqrt(A.friction_coeff * A.friction_coeff + B.friction_coeff * B.friction_coeff);
            Eigen::Vector<real, 3> f_n = calcForce(hit);
            real force = f_n.norm();
            Eigen::Vector<real, 3> n = f_n.normalized();

            Eigen::Vector<real, 3> tangent = rel_vel_new - (rel_vel_new.dot(n) * n);
            tangent.normalize();
            Eigen::Vector<real, 3> normal_force = f_n.dot(n) * n;

            real friction_mag = (rel_vel_new.dot(tangent)) / (1 / A.getMass() + 1 / B.getMass());

            Eigen::Vector<real, 3> friction_impulse = friction_mag * -tangent;

            a_force_out = friction_impulse;
            b_force_out = -friction_impulse;
            a_torque_out = calcTorque(friction_impulse, hit.A, A.future_state.getTranslation().cast<real>());
            b_torque_out = calcTorque(Eigen::Vector<real, 3>(-friction_impulse), hit.B, B.future_state.getTranslation().cast<real>());

            // if (a_force_out.hasNaN() || b_force_out.hasNaN() || a_torque_out.hasNaN() || b_torque_out.hasNaN()) {
            //     printf("Invalid friction\n");
            // }
        }

        template<typename real>
        void calcFriction(std::vector<Delta2::collision::Contact<real>> hits, const Delta2::Particle& A, const Delta2::Particle& B, Eigen::Vector<real, 3>& a_force_out, Eigen::Vector<real, 3>& a_torque_out, Eigen::Vector<real, 3>& b_force_out, Eigen::Vector<real, 3>& b_torque_out) {
            if (std::abs(A.future_state.getTime() - B.future_state.getTime()) > 1e-8) {
                throw std::invalid_argument("Input particles must have futures states with the same time");
            }
            
            for (int h = 0; h < hits.size(); h++) {
                Delta2::collision::Contact<real> hit = hits[h];

                Eigen::Vector<real, 3> a_rel_vel_new = A.current_state.pointVelocity(hit.A, A.future_state);
                Eigen::Vector<real, 3> b_rel_vel_new = B.current_state.pointVelocity(hit.B, B.future_state);

                Eigen::Vector<real, 3> rel_vel_new = a_rel_vel_new - b_rel_vel_new;

                real mu = sqrt(A.friction_coeff * A.friction_coeff + B.friction_coeff * B.friction_coeff);
                Eigen::Vector<real, 3> f_n = calcForce(hit);
                real force = f_n.norm();
                Eigen::Vector<real, 3> n = f_n.normalized();

                Eigen::Vector<real, 3> tangent = rel_vel_new - (rel_vel_new.dot(n) * n);
                tangent.normalize();
                Eigen::Vector<real, 3> normal_force = f_n.dot(n) * n;

                real friction_mag = (rel_vel_new.dot(tangent)) / (1 / A.getMass() + 1 / B.getMass());

                Eigen::Vector<real, 3> friction_impulse = friction_mag * -tangent;
                
                a_force_out += friction_impulse;
                b_force_out += -friction_impulse;
                a_torque_out += calcTorque(friction_impulse, hit.A, A.future_state.getTranslation());
                b_torque_out += calcTorque(-friction_impulse, hit.B, B.future_state.getTranslation());
            }
        }

        void friction_solve(ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size);
        void friction_solve(ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<State> future_states, std::vector<double> step_size);
    }
}

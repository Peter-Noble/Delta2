#include "forces.h"
#include "../collision_detection/contact.h"
#include "integrate.h"
#include <ittnotify.h>
#include "particle_handler.h"

using namespace Delta2;

void Delta2::model::friction_solve(ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size) {
    std::vector<Eigen::Vector3d> start_rel_vels;
    start_rel_vels.resize(hits.size());
    std::vector<double> max_force;
    max_force.resize(hits.size());

    std::vector<Eigen::Vector3d> friction_forces;
    friction_forces.resize(forces.size());
    std::vector<Eigen::Vector3d> friction_torques;
    friction_torques.resize(forces.size());

    std::vector<double> friction_factor;
    friction_factor.resize(hits.size());

    std::vector<Eigen::Vector3d> hit_friction;
    hit_friction.resize(hits.size());

    for (int i = 0; i < forces.size(); i++) {
        friction_forces[i] = {0, 0, 0};
        friction_torques[i] = {0, 0, 0};
    }
    for (int i = 0; i < hits.size(); i++) {
        friction_factor[i] = 0.1;
        hit_friction[i] = {0.0, 0.0, 0.0};
    }

    double adjust_by = 0.25;
    double rate = 0.9;
    
    int iterations = 0;
    while (iterations < 100 && hits.size() > 0) {
        for (Particle* p : particles) {
            if (!p->is_static) {
                Eigen::Vector3d force = {0.0, 0.0, 0.0};
                Eigen::Vector3d torque = {0.0, 0.0, 0.0};
                uint32_t id = particles.getLocalID(p->id);
                if (counts[id]) {
                    force = (forces[id] + friction_forces[id]) / double(counts[id]);
                    torque = (torques[id] + friction_torques[id]) / double(counts[id]);
                }
                Delta2::model::integrateEuler(*p, force + external_force(*p), torque, step_size, false);
            }
        }

        for (int i = 0; i < forces.size(); i++) {
            friction_forces[i] = {0, 0, 0};
            friction_torques[i] = {0, 0, 0};
        }
        
        for (int h_i = 0; h_i < hits.size(); h_i++) {
            collision::Contact<double> c = hits[h_i];
            Delta2::Particle& A = *c.p_a;
            Delta2::Particle& B = *c.p_b;

            if (iterations == 0) {
                double mu = sqrt(A.friction_coeff * A.friction_coeff + B.friction_coeff * B.friction_coeff);
                // Eigen::Vector<double, 3> f_n = Delta2::model::calcForce(c, 1.0, step_size);
                // max_force[h_i] = f_n.norm() * mu;
                max_force[h_i] = c.force_mag * mu;
            }

            Eigen::Vector3d a_hit_point = c.A.cast<double>();
            Eigen::Vector3d b_hit_point = c.B.cast<double>();

            Eigen::Vector<double, 3> a_rel_vel_new = A.futurePointVelocity(a_hit_point);
            Eigen::Vector<double, 3> b_rel_vel_new = B.futurePointVelocity(b_hit_point);
            Eigen::Vector<double, 3> rel_vel_new = a_rel_vel_new - b_rel_vel_new;

            //printf("%i: %f\n", h_i, rel_vel_new.norm());

            Eigen::Vector<double, 3> hit_normal = a_hit_point - b_hit_point;

            Eigen::Vector<double, 3> tangent = rel_vel_new - (rel_vel_new.dot(hit_normal.normalized()) * hit_normal.normalized());            

            Eigen::Vector3d hit_force_adjust = tangent.normalized() * max_force[h_i] * adjust_by;

            hit_friction[h_i] += hit_force_adjust;

            if (hit_friction[h_i].norm() > max_force[h_i]) {
                hit_friction[h_i] *= max_force[h_i] / hit_friction[h_i].norm();
            }

            uint32_t a_id_local = particles.getLocalID(c.p_a->id);
            uint32_t b_id_local = particles.getLocalID(c.p_b->id);

            friction_forces[a_id_local] += -hit_friction[h_i];
            friction_torques[a_id_local] += Delta2::model::calcTorque(Eigen::Vector3d(-hit_friction[h_i]), a_hit_point, A.future_state.getTranslation());
            friction_forces[b_id_local] += hit_friction[h_i];
            friction_torques[b_id_local] += Delta2::model::calcTorque(Eigen::Vector3d(hit_friction[h_i]), b_hit_point, B.future_state.getTranslation());
        }

        adjust_by *= rate;
        iterations++;
    }

    for (int p_i = 0; p_i < forces.size(); p_i++) {
        forces[p_i] += friction_forces[p_i];
        torques[p_i] += friction_torques[p_i];
    }
}

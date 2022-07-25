#include "friction.h"
#include "forces.h"
#include "../collision_detection/contact.h"
#include "integrate.h"

using namespace Delta2;


model::FrictionSolver::FrictionSolver(model::ParticleHandler& particles,
                                      std::vector<collision::Contact<double>>& input_hits) {
    ph = &particles;
    hits = &input_hits;

    max_force.resize(hits->size());

    friction_forces.resize(particles.size());
    friction_torques.resize(particles.size());

    friction_factor.resize(hits->size());

    hit_friction.resize(hits->size());

    for (int i = 0; i < particles.size(); i++) {
        friction_forces[i] = {0, 0, 0};
        friction_torques[i] = {0, 0, 0};
    }

    for (int i = 0; i < hits->size(); i++) {
        friction_factor[i] = 0.1;
        hit_friction[i] = {0.0, 0.0, 0.0};
    }

    adjust_by = 0.1;
    rate = 0.95;
    iterations = 0;
}

void model::FrictionSolver::solve(int add_iterations,
                                  std::vector<Eigen::Vector3d>& forces,
                                  std::vector<Eigen::Vector3d>& torques,
                                  std::vector<State>& future_states,
                                  double step_size) {
    int start_iteration = iterations;
    int up_to_iterations = iterations + add_iterations;
    while (iterations < up_to_iterations && hits->size() > 0) {
        for (Particle* p : *ph) {
            if (!p->is_static) {
                uint32_t id = ph->getLocalID(p->id);
                Eigen::Vector3d force = (forces[id] + friction_forces[id]);
                Eigen::Vector3d torque = (torques[id] + friction_torques[id]);

                State SState = (*ph)[id].updateState(step_size, force, torque);
                future_states[id] = SState;
            }
        }

        for (int i = 0; i < forces.size(); i++) {
            friction_forces[i] = {0, 0, 0};
            friction_torques[i] = {0, 0, 0};
        }
        
        for (int h_i = 0; h_i < hits->size(); h_i++) {
            collision::Contact<double> c = (*hits)[h_i];
            Delta2::Particle& A = *c.p_a;
            Delta2::Particle& B = *c.p_b;

            if (iterations == start_iteration) {
                double mu = sqrt(A.friction_coeff * A.friction_coeff + B.friction_coeff * B.friction_coeff);
                max_force[h_i] = c.force_mag * mu;
            }

            Eigen::Vector3d a_hit_point = c.A.cast<double>();
            Eigen::Vector3d b_hit_point = c.B.cast<double>();

            uint32_t a_id_local = ph->getLocalID(c.p_a->id);
            uint32_t b_id_local = ph->getLocalID(c.p_b->id);

            Eigen::Vector<double, 3> a_rel_vel_new = A.pointVelocity(a_hit_point, future_states[a_id_local]);
            Eigen::Vector<double, 3> b_rel_vel_new = B.pointVelocity(b_hit_point, future_states[b_id_local]);
            Eigen::Vector<double, 3> rel_vel_new = a_rel_vel_new - b_rel_vel_new;

            Eigen::Vector<double, 3> hit_normal = a_hit_point - b_hit_point;

            Eigen::Vector<double, 3> tangent = rel_vel_new - (rel_vel_new.dot(hit_normal.normalized()) * hit_normal.normalized());

            if (tangent.norm() > 1e-8) {
                Eigen::Vector3d hit_force_adjust = tangent.normalized() * max_force[h_i] * adjust_by;

                hit_friction[h_i] += hit_force_adjust;

                if (hit_friction[h_i].norm() > max_force[h_i]) {
                    hit_friction[h_i] *= max_force[h_i] / hit_friction[h_i].norm();
                }

                friction_forces[a_id_local] += -hit_friction[h_i];
                friction_torques[a_id_local] += Delta2::model::calcTorque(Eigen::Vector3d(-hit_friction[h_i]), a_hit_point, future_states[a_id_local].getTranslation());
                friction_forces[b_id_local] += hit_friction[h_i];
                friction_torques[b_id_local] += Delta2::model::calcTorque(Eigen::Vector3d(hit_friction[h_i]), b_hit_point, future_states[b_id_local].getTranslation());

                if (!A.is_static) {
                    uint32_t id = ph->getLocalID(A.id);
                    Eigen::Vector3d force = (forces[id] + friction_forces[id]);
                    Eigen::Vector3d torque = (torques[id] + friction_torques[id]);

                    State SState = (*ph)[id].updateState(step_size, force, torque);
                    future_states[id] = SState;
                }

                if (!B.is_static) {
                    uint32_t id = ph->getLocalID(B.id);
                    Eigen::Vector3d force = (forces[id] + friction_forces[id]);
                    Eigen::Vector3d torque = (torques[id] + friction_torques[id]);

                    State SState = (*ph)[id].updateState(step_size, force, torque);
                    future_states[id] = SState;
                }
            }
        }

        adjust_by *= rate;
        iterations++;
    }
}

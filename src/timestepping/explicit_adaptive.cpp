#include "explicit_adaptive.h"
#include "../collision_detection/broad_phase_embree.h"
#include "../collision_detection/full_tree_comparison.h"
#include "../model/forces.h"
#include "../model/integrate.h"

using namespace Delta2;

timestepping::ExplicitAdaptiveScheme::ExplicitAdaptiveScheme(std::vector<Particle>* ps, std::function<Eigen::Vector3d(const Particle&)> external_force, common::Options& o) : timestepping::TimestepScheme(ps, external_force, o) {
    _last_time_step_size = -1.0;
};

struct ExplicitState {
    ExplicitState() {
        force = {0, 0, 0};
        torque = {0, 0, 0};
    }
    Eigen::Vector3d force;
    Eigen::Vector3d torque;
};

void timestepping::ExplicitAdaptiveScheme::step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& view_draws) {
    double current_time = 0.0;

    if (_last_time_step_size <= 0.0) {
        _last_time_step_size = time;
    }

    while (current_time < time) {
        std::vector<std::tuple<int, int, collision::Contact<double>>> hits;

        std::vector<Eigen::Vector3d> forces;
        std::vector<Eigen::Vector3d> torques;
        std::vector<int> counts;

        bool repeat_broadphase = true;

        double time_step_size = std::min(_last_time_step_size * 1.1, time);

        collision::BroadPhaseCollisions B;

        for (int p_i = 0; p_i < _particles->size(); p_i++) {
            forces.emplace_back(0.0, 0.0, 0.0);
            torques.emplace_back(0.0, 0.0, 0.0);
            counts.emplace_back(0);
            (*_particles)[p_i].projectFutureState(time_step_size);
        }

        while (repeat_broadphase) {
            repeat_broadphase = false;
            
            B = collision::broadPhaseEmbree(*_particles, time_step_size);

            for (collision::BroadPhaseCollision& b : B) {
                int a_id = b.first.first; // geo id
                int b_id = b.second.first;

                Eigen::Vector3d a_centre = (*_particles)[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = (*_particles)[b_id].current_state.getTranslation();

                // common::Viewer view;
                // view.addParticleInterval((*_particles)[a_id]);
                // view.addParticleInterval((*_particles)[b_id]);
                // view.show();

                double min_time;
                std::vector<collision::ContinuousContact<double>> Cs = collision::compareTreesFullContinuous<double, 8, 8>((*_particles)[a_id], (*_particles)[b_id], time_step_size, min_time);
                double old_time_step_size = time_step_size;

                if (Cs.size() > 0) {
                    double min_scaling_seen = 1.0;
                    for (collision::ContinuousContact<double>& c : Cs) {
                        Eigen::Vector3d hit_normal = (c.A - c.B);
                        if (hit_normal.norm() < 1e-8) {
                            if (c.toc == 0.0) {
                                throw std::runtime_error("Invalid configuration");
                            }
                            repeat_broadphase = true;
                            min_scaling_seen = std::min(min_scaling_seen, c.toc * 0.75);
                            continue;
                        }
                        double interaction_dist = c.eps_a + c.eps_b;

                        Eigen::Vector<double, 3> a_vel = c.p_a->futurePointVelocity(c.A);
                        Eigen::Vector<double, 3> b_vel = c.p_b->futurePointVelocity(c.B);

                        Eigen::Vector<double, 3> rel_vel = a_vel - b_vel;

                        if (hit_normal.normalized().dot(rel_vel.normalized()) <= 0.0) {
                            if (hit_normal.norm() <= interaction_dist) {
                                // End point is inside range (since it's moving towards each other) => compute start point
                                double proj_vel = -(hit_normal.dot(rel_vel)/hit_normal.dot(hit_normal)) * hit_normal.norm(); // +ve
                                double depth = interaction_dist - hit_normal.norm();
                                double proj_dist = proj_vel * c.toc * old_time_step_size;
                                if (depth < proj_dist) {
                                    // This is the first frame penetrating this eps boundry
                                    min_scaling_seen = std::min(min_scaling_seen, c.toc * (2.0 - depth / proj_dist) / 2.0);
                                }
                                else {
                                    // The contact point was already inside the eps boundry at the start of the timestep
                                    if (hit_normal.norm() / proj_dist < 4.0) {
                                        min_scaling_seen = std::min(min_scaling_seen, 0.5 * c.toc);
                                    }
                                }
                            }
                            // Else there is not point in the range that interacts => don't make the timestep smaller.

                            // min_scaling_seen = std::min(min_scaling_seen, c.toc * (hit_normal.norm() / (2.0 * interaction_dist) + 0.5));
                        }
                    }
                    time_step_size *= min_scaling_seen;
                    for (int p_i = 0; p_i < _particles->size(); p_i++) {
                        (*_particles)[p_i].projectFutureState(time_step_size);
                    }

                }
            }
        }

        _last_time_step_size = time_step_size;

        for (Particle& p : *_particles) {
            p.future_state = p.current_state;
        }

        for (collision::BroadPhaseCollision& b : B) {
            int a_id = b.first.first; // geo id
            int b_id = b.second.first;

            Eigen::Vector3d a_centre = (*_particles)[a_id].current_state.getTranslation();
            Eigen::Vector3d b_centre = (*_particles)[b_id].current_state.getTranslation();

            std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>((*_particles)[a_id], (*_particles)[b_id]);
            std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

            for (collision::Contact<double>& c : filtered) {
                hits.push_back(std::make_tuple(a_id, b_id, c));
            }
        }
        
        for (auto& [a_id, b_id, c] : hits) {
            Eigen::Vector3d a_centre = (*_particles)[a_id].current_state.getTranslation();
            Eigen::Vector3d b_centre = (*_particles)[b_id].current_state.getTranslation();

            Eigen::Vector3d a_hit_point = c.A.cast<double>();
            Eigen::Vector3d b_hit_point = c.B.cast<double>();

            double mass_multiply = 0.0;
            if (!(*_particles)[a_id].is_static) {
                mass_multiply += (*_particles)[a_id].getMass();
            }
            if ((*_particles)[b_id].is_static) {
                mass_multiply += (*_particles)[b_id].getMass();
            }
            Eigen::Vector3d f = model::calcForce(c).cast<double>() * mass_multiply;
            Eigen::Vector3d t_a = model::calcTorque(f, a_hit_point, a_centre);
            Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), b_hit_point, b_centre);

            forces[a_id] += f;
            torques[a_id] += t_a;
            counts[a_id]++;
            forces[b_id] += -f;
            torques[b_id] += t_b;
            counts[b_id]++;

            view_draws.push_back(std::make_pair(a_hit_point, b_hit_point));
        }

        for (int p_i = 0; p_i < _particles->size(); p_i++) {
            Particle& p = (*_particles)[p_i];
            if (!p.is_static) {
                Eigen::Vector3d force = {0.0, 0.0, 0.0};
                Eigen::Vector3d torque = {0.0, 0.0, 0.0};
                if (counts[p_i]) {
                    force = forces[p_i] / double(counts[p_i]);
                    torque = torques[p_i] / double(counts[p_i]);
                }
                model::integrateEuler(p, force + _external_force(p), torque, _last_time_step_size, false);
            } else {
                model::staticApply(p, _last_time_step_size, false);
            }
        }

        // printf("Pre friction force:  {%f, %f, %f}\n", forces[0].x(), forces[0].y(), forces[0].z());
        // printf("Pre friction torque: {%f, %f, %f}\n", torques[0].x(), torques[0].y(), torques[0].z());
        // for (auto& [a_id, b_id, c] : hits) {
        //     Eigen::Vector3d a_hit_point = c.A.cast<double>();
        //     Eigen::Vector3d b_hit_point = c.B.cast<double>();

        //     Eigen::Vector3d friction_f_a, friction_f_b, friction_t_a, friction_t_b;
        //     model::calcFriction(c, (*_particles)[a_id], (*_particles)[b_id], friction_f_a, friction_t_a, friction_f_b, friction_t_b);

        //     view_draws.emplace_back(c.A, c.A + friction_f_a.normalized());

        //     torques[b_id] += friction_t_b;
        //     torques[a_id] += friction_t_a;
        //     forces[b_id] += friction_f_b;
        //     forces[a_id] += friction_f_a;
        // }
        Delta2::model::friction_solve(*_particles, hits, forces, torques, counts, _external_force, _last_time_step_size);

        // printf("Post friction force:  {%f, %f, %f}\n", forces[0].x(), forces[0].y(), forces[0].z());
        // printf("Post friction torque: {%f, %f, %f}\n", torques[0].x(), torques[0].y(), torques[0].z());

        for (int p_i = 0; p_i < _particles->size(); p_i++) {
            Particle& p = (*_particles)[p_i];
            if (!p.is_static) {
                Eigen::Vector3d force = {0.0, 0.0, 0.0};
                Eigen::Vector3d torque = {0.0, 0.0, 0.0};
                if (counts[p_i]) {
                    force = forces[p_i] / double(counts[p_i]);
                    torque = torques[p_i] / double(counts[p_i]);
                }
                model::integrateEuler(p, force + _external_force(p), torque, _last_time_step_size);
            } else {
                model::staticApply(p, _last_time_step_size);
            }
        }
        current_time += _last_time_step_size;
        printf("Step to: %f\n", current_time);
    }
}

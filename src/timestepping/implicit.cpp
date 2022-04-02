#include "implicit.h"
#include "../collision_detection/broad_phase_embree.h"
#include "../collision_detection/full_tree_comparison.h"
#include "../model/forces.h"
#include "../model/integrate.h"

using namespace Delta2;

timestepping::ImplicitScheme::ImplicitScheme(std::vector<Particle>* ps, std::function<Eigen::Vector3d(const Particle&)> external_force, common::Options& o) : timestepping::TimestepScheme(ps, external_force, o) {};

void timestepping::ImplicitScheme::step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& view_draws) {
    for (Particle& p : *_particles) {
        p.future_state = p.current_state;
    }

    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<int> counts;
    for (int p_i = 0; p_i < _particles->size(); p_i++) {
        forces.emplace_back(0.0, 0.0, 0.0);
        torques.emplace_back(0.0, 0.0, 0.0);
        counts.emplace_back(0);
    }

    collision::BroadPhaseCollisions B = collision::broadPhaseEmbree(*_particles, time);

    std::vector<std::tuple<int, int, collision::Contact<double>>> hits;

    int num_particles = _particles->size();
    Eigen::SparseMatrix<int> interaction_graph(num_particles, num_particles);

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
            model::integrateEuler(p, force + _external_force(p) , torque, time, false);
        } else {
            model::staticApply(p, time, false);
        }
    }

    for (auto& [a_id, b_id, c] : hits) {
        Eigen::Vector3d a_hit_point = c.A.cast<double>();
        Eigen::Vector3d b_hit_point = c.B.cast<double>();

        Eigen::Vector3d friction_f_a, friction_f_b, friction_t_a, friction_t_b;
        model::calcFriction(c, (*_particles)[a_id], (*_particles)[b_id], friction_f_a, friction_t_a, friction_f_b, friction_t_b);

        torques[b_id] += 100.0 * friction_t_b;
        torques[a_id] += 100.0 * friction_t_a;
        forces[b_id] += 100.0 * friction_f_b;
        forces[a_id] += 100.0 * friction_f_a;
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
            model::integrateEuler(p, force + _external_force(p), torque, time);
        } else {
            model::staticApply(p, time);
        }
    }
}

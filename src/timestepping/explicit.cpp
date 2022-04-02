#include "explicit.h"
#include "../collision_detection/broad_phase_embree.h"
#include "../collision_detection/full_tree_comparison.h"
#include "../model/forces.h"
#include "../model/integrate.h"

using namespace Delta2;

timestepping::ExplicitScheme::ExplicitScheme(std::vector<Particle>* ps, common::Options& o) : timestepping::TimestepScheme(ps, o) {
    _step_size = 0.01;
    for (int i = 0; i < ps->size(); i++) {
        _global_to_local_ids[i] = i;
    }
};

void timestepping::ExplicitScheme::stepSetup() {
    for (Particle& p : *_particles) {
        p.future_state = p.current_state;
    }
    _forces.clear();
    _torques.clear();
    _counts.clear();
    for (int p_i = 0; p_i < _particles->size(); p_i++) {
        _forces.emplace_back(0.0, 0.0, 0.0);
        _torques.emplace_back(0.0, 0.0, 0.0);
        _counts.emplace_back(0);
    }
}

void timestepping::ExplicitScheme::broadPhase() {
    _broad_phase_collisions = collision::broadPhaseEmbree(*_particles, _step_size);
}

void timestepping::ExplicitScheme::narrowPhase(collision::BroadPhaseCollision& b) {
    int a_id = b.first.first; // geo id
    int b_id = b.second.first;
    
    std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>((*_particles)[a_id], (*_particles)[b_id]);
    std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

    for (collision::Contact<double>& c : filtered) {
        _hits.push_back(std::make_tuple(a_id, b_id, c));
    }
}


void timestepping::ExplicitScheme::solveContacts() {
    for (auto& [a_id, b_id, c] : _hits) {
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

        _forces[a_id] += f;
        _torques[a_id] += t_a;
        _counts[a_id]++;
        _forces[b_id] += -f;
        _torques[b_id] += t_b;
        _counts[b_id]++;

        _view_draws.push_back(std::make_pair(a_hit_point, b_hit_point));
    }
}

void timestepping::ExplicitScheme::integrate(u_int32_t p_i, bool apply) {
    Particle& p = (*_particles)[p_i];
    if (!p.is_static) {
        Eigen::Vector3d force = {0.0, 0.0, 0.0};
        Eigen::Vector3d torque = {0.0, 0.0, 0.0};
        if (_counts[p_i] > 0) {
            force = _forces[p_i] / _counts[p_i];
            torque = _torques[p_i] / _counts[p_i];
        }
        model::integrateEuler(p, force + external_force(p), torque, _step_size, false);
    } else {
        model::staticApply(p, _step_size, false);
    }
}

void timestepping::ExplicitScheme::friction() {
    Delta2::model::friction_solve(*_particles, _hits, _forces, _torques, _counts, external_force, _step_size, _global_to_local_ids);
}

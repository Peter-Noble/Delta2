#include "timestep_scheme.h"

using namespace Delta2;

timestepping::TimestepScheme::TimestepScheme(std::vector<Particle>* p, common::Options& o) : _particles(p), _opt(o) {
    for (int p_i = 0; p_i < _particles->size(); p_i++) {
        (*_particles)[p_i].id = p_i;
    }
    external_force = default_gravity_callback;
};

timestepping::TimestepScheme::~TimestepScheme() {};

void timestepping::TimestepScheme::step() {
    _view_draws.clear();

    stepSetup();
    broadPhase();
    narrowPhase();

    for (int p_i = 0; p_i < _particles->size(); p_i++) {
        integrate(p_i, false);
    }

    friction();

    for (int p_i = 0; p_i < _particles->size(); p_i++) {
        integrate(p_i, true);
    }
}

void timestepping::TimestepScheme::setStepSize(double step) {
    if (step <= 0.0) {
        throw std::invalid_argument("Step size must be positive");
    }
    _step_size = step;
}

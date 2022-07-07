#pragma once

#include <Eigen/Dense>

#include "../../common/cli_options.h"
#include "../../model/particle_handler.h"
#include "../../collision_detection/contact.h"

namespace Delta2 {
    namespace strategy {
        class FrictionStrategy {
        public:
            FrictionStrategy(common::Options& opt);
            virtual void solve(model::ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size) = 0;
            virtual void solve(model::ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<State> future_states, std::vector<double> step_size) = 0;
        protected:
            common::Options _opt;
        };
    }
}
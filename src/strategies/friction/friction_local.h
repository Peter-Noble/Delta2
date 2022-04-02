#pragma once

#include "friction_strategy.h"

namespace Delta2 {
    namespace strategy {
        class FrictionLocal : public FrictionStrategy {
        public:
            FrictionLocal(common::Options& opt);
            void solve(model::ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size) override;
        };
    }
}
#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void start_intersecting(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation({0, 0, 0});
            }

            {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation({1.95, 1.95, 0});
            }
        }
    }
}

#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"

namespace Delta2 {
    namespace scenarios {
        void pairs(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            force_strategy.external_force = [](const Delta2::Particle& p) {
                return Eigen::Vector3d({0, 0, 0});
            };

            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Delta2::common::sphere(1, 6, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            int scenario_size = 10;
            if (globals::opt.scenario_size >= 0) {
                scenario_size = globals::opt.scenario_size;
            }

            const double offset = 5;

            for (int x = 0; x < scenario_size; x++) {
                double separation = globals::opt.rand_float(1.9) + 0.1;
                {
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation(x * Eigen::Vector3d({offset, 0, 0}) + Eigen::Vector3d({0.5, -2, 0}));
                    p.current_state.setVelocity({0, separation, 0});
                }
                {
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation(x * Eigen::Vector3d({offset, 0, 0}) + Eigen::Vector3d({-0.5, 2, 0}));
                    p.current_state.setVelocity({0, -separation, 0});
                }
            }
        }       
    }
}

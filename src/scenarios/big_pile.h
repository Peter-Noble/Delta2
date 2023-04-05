#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void pyramid(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            
            Delta2::globals::opt.seed_rand(12345, Delta2::globals::opt.scenario_size);
            {
                Delta2::common::plane(100, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.25);
                p.is_static = true;
            }

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));
            
            int size = globals::opt.scenario_size;

            for (int z = 0; z < size; z++) {
                for (int x = 0; x < size-z; x++) {
                    for (int y = 0; y < size-z; y++) {
                        auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
                        double offset = -((double)(size-1-z))/2.0;
                        p.current_state.setTranslation({(offset + x) * 2.01, (offset + y) * 2.01, 1.01 + 2.1 * z + Delta2::globals::opt.rand_float(0.1)});
                    }
                }
            }
        }
    }
}

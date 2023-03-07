#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void fast_moving(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            int scenario_size = 1;
            if (globals::opt.scenario_size >= 0) {
                scenario_size = globals::opt.scenario_size;
            }
            
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            
            {
                Delta2::common::plane(20, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.is_static = true;
            }

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            for (int i = 0; i < scenario_size; i++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation({globals::opt.rand_float(10.0)-5, globals::opt.rand_float(10.0)-5, 10 + 3 * i});
                p.current_state.setVelocity({0, 0, -50});
                p.current_state.setAngular({globals::opt.rand_float(10.0), globals::opt.rand_float(10.0), globals::opt.rand_float(10.0)});
            }
        }
    }
}

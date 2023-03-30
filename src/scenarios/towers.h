#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void towers(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            int max_x = 1;
            if (globals::opt.scenario_size >= 0) {
                max_x = globals::opt.scenario_size;
            }
            const double spacing_x = 8.0;
            const int max_y = 4;
            const double spacing_y = 30.0;
            const int max_i = 20;
            const int tilt_offset = 3;
            
            {
                Delta2::common::plane(std::max(10.0, std::max(max_x * spacing_x, max_y * spacing_y * 2)), V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.is_static = true;
            }

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            for (int y = 0; y < max_y; y++) {
                for (int x = 0; x < max_x; x++) {
                    for (int i = 0; i < max_i; i++) {
                        auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                        p.current_state.setTranslation({-(max_x / 2 * spacing_x) + x * spacing_x, -(max_y / 2 * spacing_y) + y * spacing_y + (i * -y * tilt_offset) * 0.2, 1.1 + i * 2.1});
                        p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({0.0, 0.0, (double)globals::opt.rand_float(igl::PI*2.0)})));
                    }
                }
            }
        }
    }
}

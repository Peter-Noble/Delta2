#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void round_tower(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Delta2::common::plane(100, V, F);
            std::shared_ptr<Delta2::MeshData> P(new Delta2::MeshData(V, F, globals::opt, true));

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            std::shared_ptr<Delta2::MeshData> S(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/deformed_icosphere.obj", globals::opt));

            int max_z = 1;
            if (globals::opt.scenario_size >= 0) {
                max_z = globals::opt.scenario_size;
            }

            int size = 1;
            if (globals::opt.scenario_size2 >= 0) {
                size = globals::opt.scenario_size2;
            }

            const double r = 7*2;
            const int ring_num = 18*2;

            for (int i = 0; i < size; i++) {
                Eigen::Vector3d offset({100 * i, 0, 0});

                {
                    auto& p = particles.emplace_back(P, 1.0, 1.0, globals::opt.geo_eps);
                    p.is_static = true;
                    p.current_state.setTranslation(offset);
                }

                for (int z = 0; z < max_z; z++) {
                    for (int x = 0; x < ring_num; x++) {
                        const Eigen::Quaterniond rot = common::eulerAnglesToQuaternion(Eigen::Vector3d({0, 0, common::degToRad((x+z/2.0)*360.0/ring_num)}));
                        Eigen::Vector3d pos_r = common::transform(Eigen::Vector3d({r, 0, 0}), rot);
                        auto& p = particles.emplace_back(M, 4.0, 1.0, globals::opt.geo_eps);  // heavy cube (Aluminium) 2700
                        p.current_state.setTranslation(offset + pos_r + Eigen::Vector3d({0, 0, z * 2.05 + 1.025}));
                        p.current_state.setRotation(rot);
                    }
                }

                {
                    auto& p = particles.emplace_back(S, 50.0, 1.0, globals::opt.geo_eps);
                    p.current_state.setTranslation(offset + Eigen::Vector3d({r + 10, 0, 6}));
                    p.current_state.setVelocity({-100, 0, 3});
                }
            }
        }
    }
}

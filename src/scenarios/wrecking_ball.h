#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void wrecking_ball(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            
            {
                Delta2::common::plane(100, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.25);
                p.is_static = true;
            }

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            int max_z = 1;
            if (globals::opt.scenario_size >= 0) {
                max_z = globals::opt.scenario_size;
            }

            const double r = 7;
            const int ring_num = 18;
            for (int z = 0; z < max_z; z++) {
                for (int x = 0; x < ring_num; x++) {
                    const Eigen::Quaterniond rot = common::eulerAnglesToQuaternion(Eigen::Vector3d({0, 0, common::degToRad((x+z/2.0)*360.0/ring_num)}));
                    Eigen::Vector3d pos_r = common::transform(Eigen::Vector3d({r, 0, 0}), rot);
                    auto& p = particles.emplace_back(M, 1.0, 1.0, 0.25);  // heavy cube (Aluminium) 2700
                    p.current_state.setTranslation(pos_r + Eigen::Vector3d({0, 0, z * 2.05 + 1.025}));
                    p.current_state.setRotation(rot);
                }
            }

            std::shared_ptr<Delta2::MeshData> C(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/chainlink.obj", globals::opt));
            
            {
                auto& p = particles.emplace_back(C, 1.0, 1.0, 0.25);
                p.current_state.setTranslation({0, 0, max_z * 2.05 + 4 + 10});
                p.is_static = true;
            }

            const int chain_length = 10;

            for (int l = 0; l < chain_length; l++) {
                const Eigen::Quaterniond rot = common::eulerAnglesToQuaternion(Eigen::Vector3d({common::degToRad(90.0) * (l + 1), common::degToRad(90.0), 0}));
                auto& p = particles.emplace_back(C, 50.0, 1.0, 0.25);
                p.current_state.setTranslation({1 + l * 3.0, 0, max_z * 2.05 + 4 + 10 - 2});
                p.current_state.setRotation(rot);
            }

            // std::shared_ptr<Delta2::MeshData> W(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/wreckingball.obj", globals::opt));
            
            // {
            //     const Eigen::Quaterniond rot = common::eulerAnglesToQuaternion(Eigen::Vector3d({common::degToRad(90.0) * (chain_length + 1), common::degToRad(-90.0), common::degToRad(0)}));
            //     auto& p = particles.emplace_back(W, 1.0, 1.0, 0.25);
            //     p.current_state.setTranslation({(chain_length + 1) * 3.0, 0, max_z * 2.05 + 4 + 10 - 2});
            //     p.current_state.setRotation(rot);
            // }
        }
    }
}

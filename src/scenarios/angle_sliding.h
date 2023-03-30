#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void angle_sliding(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            
            Delta2::globals::opt.seed_rand(12345, Delta2::globals::opt.scenario_size);
            {
                Delta2::common::plane(100, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.25);
                p.is_static = true;
            }

            {
                Delta2::common::cube(1.0, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
                p.current_state.setVelocity({Delta2::globals::opt.rand_float(3), Delta2::globals::opt.rand_float(0.5), Delta2::globals::opt.rand_float(0.5)});
                p.current_state.setTranslation({0, 0, 1.3 + Delta2::globals::opt.rand_float(0.5)});
                p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({Delta2::common::degToRad(180 - Delta2::globals::opt.rand_float(360)), Delta2::common::degToRad(180 - Delta2::globals::opt.rand_float(360)), Delta2::common::degToRad(180 - Delta2::globals::opt.rand_float(360))})));
                p.current_state.setAngular({5-Delta2::globals::opt.rand_float(10), 5-Delta2::globals::opt.rand_float(10), 5-Delta2::globals::opt.rand_float(10)});
            }
        }

        void angle_sliding_first(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            
            Delta2::globals::opt.seed_rand(12345, Delta2::globals::opt.scenario_size);
            {
                Delta2::common::plane(100, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.25);
                p.is_static = true;
            }

            {
                Delta2::common::cube(1.0, V, F);
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));
                auto& p = particles.emplace_back(M, 1.0, 1.0, 0.1);
                p.current_state.setRotation({-0.060924, -0.52497188, -0.8372, 0.1404});
                p.current_state.setAngular({10.04053, -7.81168, 0.3027});
                p.current_state.setTranslation({-0.99927, -1.69319, 1.41747});
                p.current_state.setVelocity({-1.9568, -2.609366, -2.43002});
            }
        }
    }
}

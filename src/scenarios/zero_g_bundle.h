#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"

namespace Delta2 {
    namespace scenarios {
        void zero_g_bundle_orig(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            force_strategy.external_force = [](const Delta2::Particle& p) {
                return Eigen::Vector3d({0, 0, 0});
            };

            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::cube(0.5, V, F);
            std::shared_ptr<Delta2::MeshData> M1(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::sphere(4, 6, V, F);
            std::shared_ptr<Delta2::MeshData> M2(new Delta2::MeshData(V, F, globals::opt));

            const double start_offset = 10;
            const double offset = 8;
            const int line = 10;
            const double vel_mult = 10;

            Eigen::Vector3d normal = {0, 1, 0};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((10 + start_offset + x * offset) * normal + Eigen::Vector3d({0, 0, -10}));
                p.current_state.setVelocity(-normal * vel_mult);
            }

            normal = {-1, 0, 0};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation(((start_offset + 2.0) + x * offset) * normal + Eigen::Vector3d({-10, -10, -10}));
                p.current_state.setVelocity(-normal * vel_mult * 1.5);
            }

            normal = {0, 1, 1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((1.0 + start_offset + x * offset) * normal);
                p.current_state.setVelocity(-normal * vel_mult);
            }

            normal = {0, -1, 1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((2.0 + start_offset + x * offset) * normal);
                p.current_state.setVelocity(-normal * vel_mult);
            }

            normal = {0, -1, 1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((10.0 + start_offset + x * offset) * normal + Eigen::Vector3d({15, 0, 0}));
                p.current_state.setVelocity(-normal * vel_mult * 0.75);
            }

            normal = {-1, -1, -1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M2, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((10.0 + start_offset + x * offset * 2.5) * normal + Eigen::Vector3d({10, 0, 0}));
                p.current_state.setVelocity(-normal * vel_mult * 0.6);
            }

            normal = {1, 1, 1};
            normal.normalize();

            for (int i = -2; i <= 2; i++) {
                for (int j = -2; j <= 2; j++) {
                    for (int x = 0; x < line * 2; x++) {
                        auto& p = particles.emplace_back(M1, 1.0, 10.0, 0.25);
                        p.current_state.setTranslation((4.0 + start_offset + x * offset) * normal + Eigen::Vector3d({6, i * 2, j * 2}));
                        p.current_state.setVelocity(-normal * vel_mult * 2.0);
                    }
                }   
            }

            normal = {1, 0, 1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M1, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((5.0 + start_offset + x * offset) * normal + Eigen::Vector3d({10, 0, 0}));
                p.current_state.setVelocity(-normal * vel_mult * 2.0);
            }

            normal = {0.5, 1, 0};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((3.0 + start_offset + x * offset) * normal);
                p.current_state.setVelocity(-normal * vel_mult);
            }

            normal = {0, 1, -1};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((4.0 + start_offset + x * offset) * normal);
                p.current_state.setVelocity(-normal * vel_mult);
            }

            normal = {1, 0, 0};
            normal.normalize();

            for (int x = 0; x < line; x++) {
                auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                p.current_state.setTranslation((5.0 + start_offset + x * offset) * normal);
                p.current_state.setVelocity(-normal * vel_mult);
            }
        }

        void zero_g_bundle(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            force_strategy.external_force = [](const Delta2::Particle& p) {
                return Eigen::Vector3d({0, 0, 0});
            };

            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Delta2::common::cube(1.0, V, F);
            std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::cube(0.5, V, F);
            std::shared_ptr<Delta2::MeshData> M1(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::sphere(4, 6, V, F);
            std::shared_ptr<Delta2::MeshData> M2(new Delta2::MeshData(V, F, globals::opt));

            int scenario_size = 1;
            if (globals::opt.scenario_size >= 0) {
                scenario_size = globals::opt.scenario_size;
            }

            const double offset = 26;

            for (int x = 0; x < scenario_size; x++) {
                double separation = globals::opt.rand_float(1.0);

                {
                    auto& p = particles.emplace_back(M2, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset, 0, 0});
                    p.current_state.setVelocity({separation, 0, 0});
                }

                for (int i = 0; i < 4; i++) {
                    for (int y = 0; y < 4; y++) {
                        for (int z = 0; z < 4; z++) {
                            separation = globals::opt.rand_float(2.0);  
                            auto& p = particles.emplace_back(M1, 1.0, 10.0, 0.25);
                            p.current_state.setTranslation({x * offset + 5 + 5 * i, -4 + y * 1.5, -4 + z * 1.5});
                            p.current_state.setVelocity({-separation, 0, 0});
                        }
                    }
                }

                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, 7, 0});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, -7, 0});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, 0, 7});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, 0, -7});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, +5, +5});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, +5, -5});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, -5, +5});
                    p.current_state.setVelocity({separation, 0, 0});
                }
                {
                    separation = globals::opt.rand_float(4.0) - 2.0;
                    auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
                    p.current_state.setTranslation({x * offset + 2, -5, -5});
                    p.current_state.setVelocity({separation, 0, 0});
                }
            }
        }       
    }
}

#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void hoppers(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Delta2::common::plane(26, V, F);
            std::shared_ptr<Delta2::MeshData> P(new Delta2::MeshData(V, F, globals::opt, true));
            std::shared_ptr<Delta2::MeshData> H(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/hopper.obj", globals::opt));
            

            std::vector<std::shared_ptr<Delta2::MeshData>> large;
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/1.obj", globals::opt));
                large.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/2.obj", globals::opt));
                large.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/3.obj", globals::opt));
                large.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/4.obj", globals::opt));
                large.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/5.obj", globals::opt));
                large.push_back(M);
            };

            std::vector<std::shared_ptr<Delta2::MeshData>> medium;
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/med_1.obj", globals::opt));
                medium.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/med_2.obj", globals::opt));
                medium.push_back(M);
            };

            std::vector<std::shared_ptr<Delta2::MeshData>> small;
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/small_1.obj", globals::opt));
                small.push_back(M);
            };
            {
                std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(std::string(std::getenv("MESH_DIR")) + "/rocks/small_2.obj", globals::opt));
                small.push_back(M);
            };

            int size = 1;
            if (globals::opt.scenario_size >= 0) {
                size = globals::opt.scenario_size;
            }

            Delta2::common::cube(0.1, V, F);
            std::shared_ptr<Delta2::MeshData> C(new Delta2::MeshData(V, F, globals::opt));

            for (int x = 0; x < size; x++) {
                for (int y = 0; y < 1; y++) {
                    Eigen::Vector3d centre = {24.6 * x, 24.6 * y, 0};
                    {
                        auto& p = particles.emplace_back(P, 1.0, 1.0, globals::opt.geo_eps);
                        p.current_state.setTranslation(centre + Eigen::Vector3d({0, 0, -25}));
                        p.is_static = true;
                    }
                    {
                        auto& p = particles.emplace_back(H, 1.0, 1.0, globals::opt.geo_eps);
                        p.current_state.setRotation(common::eulerAnglesToQuaternion(Eigen::Vector3d({common::degToRad(90.0), 0, 0})));
                        p.current_state.setTranslation(centre);
                        p.is_static = true;
                    }
                    for (int px = -2; px < 3; px++) {
                        for (int py = -2; py < 3; py++) {
                            for (int pz = 0; pz < 2; pz++) {
                                if ((px + py + pz) % 2 != 0) {
                                    auto& p = particles.emplace_back(large[Delta2::globals::opt.rand(large.size())], 1.0, 1.0, globals::opt.geo_eps);
                                    // auto& p = particles.emplace_back(C, 1.0, 1.0, globals::opt.geo_eps);
                                    p.current_state.setTranslation(centre + Eigen::Vector3d({px * 3.0, py * 3.0, pz * 2.5}));
                                }
                            }
                            for (int pz = 0; pz < 3; pz++) {
                                if ((px + py + pz) % 2 != 0) {
                                    auto& p = particles.emplace_back(medium[Delta2::globals::opt.rand(medium.size())], 1.0, 1.0, globals::opt.geo_eps);
                                    // auto& p = particles.emplace_back(C, 1.0, 1.0, globals::opt.geo_eps);
                                    p.current_state.setTranslation(centre + Eigen::Vector3d({px * 5, py * 5, 4 + pz * 1.75}));
                                }
                            }
                            for (int pz = 0; pz < 6; pz++) {
                                if ((px + py + pz) % 2 != 0) {
                                    auto& p = particles.emplace_back(small[Delta2::globals::opt.rand(small.size())], 1.0, 1.0, globals::opt.geo_eps);
                                    // auto& p = particles.emplace_back(C, 1.0, 1.0, globals::opt.geo_eps);
                                    p.current_state.setTranslation(centre + Eigen::Vector3d({px * 4.0, py * 4.0, 10 + pz * 1}));
                                }
                            }
                        }
                    }

                    // auto& p = particles.emplace_back(C, 1.0, 1.0, 0.25);
                    // p.current_state.setTranslation(centre + Eigen::Vector3d({6.0, 4.0, -2}));
                    // p.current_state.setVelocity(Eigen::Vector3d({0, 0, -10}));
                }
            }
        }
    }
}

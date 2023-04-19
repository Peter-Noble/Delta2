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
            Delta2::common::plane(100, V, F);
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

            for (int x = 0; x < size; x++) {
                for (int y = 0; y < size; y++) {
                    Eigen::Vector3d centre = {24.6 * x, 24.6 * y, 0};
                    {
                        auto& p = particles.emplace_back(P, 1.0, 1.0, 0.25);
                        p.current_state.setTranslation(centre + Eigen::Vector3d({0, 0, -5}));
                        p.is_static = true;
                    }
                    {
                        auto& p = particles.emplace_back(H, 1.0, 1.0, 0.25);
                        p.current_state.setTranslation(centre);
                        p.is_static = true;
                    }
                    for (int px = 0; px < 3; px++) {
                        for (int py = 0; py < 3; py++) {
                            for (int pz = 0; pz < 2; pz++) {
                                auto& p = particles.emplace_back(large[Delta2::globals::opt.rand(large.size())], 1.0, 1.0, 0.25);
                                p.current_state.setTranslation()
                            }
                            for (int pz = 0; pz < 4; pz++) {
                                
                            }
                            for (int pz = 0; pz < 6; pz++) {
                                
                            }
                        }
                    }
                }
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
    }
}

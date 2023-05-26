#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "igl/PI.h"

namespace Delta2 {
    namespace scenarios {
        void waterfall(std::vector<Delta2::Particle>& particles, strategy::ContactForceStrategy& force_strategy) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            Delta2::common::plane(4, V, F);
            std::shared_ptr<Delta2::MeshData> P(new Delta2::MeshData(V, F, globals::opt, true));

            Delta2::common::sphere(1, 6, V, F);
            std::shared_ptr<Delta2::MeshData> S(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::cube(0.5, V, F);
            std::shared_ptr<Delta2::MeshData> C(new Delta2::MeshData(V, F, globals::opt));

            Delta2::common::sphere(4, 6, V, F);
            std::shared_ptr<Delta2::MeshData> S2(new Delta2::MeshData(V, F, globals::opt));

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

            int scenario_size = 1;
            if (globals::opt.scenario_size >= 0) {
                scenario_size = globals::opt.scenario_size;
            }
            bool use_hopper_particles = false;
            if (globals::opt.scenario_lod > 0) {
                use_hopper_particles = true;
            }

            const double x_spacing = 4.1;
            double x_offset = scenario_size / 2.0 * x_spacing;

            int height = 10;
            double y_z_spacing = 5.0;
            double y_z_offset = height / 2.0 * y_z_spacing;

            for (int i = 0; i < height; i++) {
                for (int x = 0; x < scenario_size; x++) {
                    {
                        auto& p = particles.emplace_back(P, 1.0, 10.0, globals::opt.geo_eps);
                        p.current_state.setTranslation(x * Eigen::Vector3d({x_spacing, 0, 0}) + Eigen::Vector3d({-x_offset, 2 - y_z_offset + i * y_z_spacing, 0.41 - y_z_offset + i * y_z_spacing}));
                        p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({2.0/360.0*igl::PI*2.0, 0, 0})));
                        p.is_static = true;
                    }
                    {
                        auto& p = particles.emplace_back(P, 1.0, 10.0, globals::opt.geo_eps);
                        p.current_state.setTranslation(x * Eigen::Vector3d({x_spacing, 0, 0}) + Eigen::Vector3d({-x_offset, -y_z_offset + 4.5 + i * y_z_spacing, -y_z_offset + 2.81 + i * y_z_spacing}));
                        p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({80.0/360.0*igl::PI*2.0, 0, 0})));
                        p.is_static = true;
                    }
                }
            }

            if (use_hopper_particles) {
                for (int px = 0; px < scenario_size; px++) {
                    for (int py = 0; py < 5; py++) {
                        for (int pz = 0; pz < 2; pz++) {
                            if ((px + py + pz) % 2 != 0) {
                                auto& p = particles.emplace_back(large[Delta2::globals::opt.rand(large.size())], 1.0, 1.0, globals::opt.geo_eps);
                                p.current_state.setTranslation(Eigen::Vector3d({px * x_spacing - x_offset, py * 3.0, pz * 2.5 + 25}));
                            }
                        }
                        for (int pz = 0; pz < 3; pz++) {
                            if ((px + py + pz) % 2 != 0) {
                                auto& p = particles.emplace_back(medium[Delta2::globals::opt.rand(medium.size())], 1.0, 1.0, globals::opt.geo_eps);
                                p.current_state.setTranslation(Eigen::Vector3d({px * x_spacing - x_offset, py * 5, 4 + pz * 1.75 + 25}));
                            }
                        }
                        for (int pz = 0; pz < 6; pz++) {
                            if ((px + py + pz) % 2 != 0) {
                                auto& p = particles.emplace_back(small[Delta2::globals::opt.rand(small.size())], 1.0, 1.0, globals::opt.geo_eps);
                                p.current_state.setTranslation(Eigen::Vector3d({px * x_spacing - x_offset, py * 4.0, 10 + pz * 1 + 25}));
                            }
                        }
                    }
                }
            }
            else {
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < scenario_size; x++) {
                        if (globals::opt.rand_float(1.0) < 0.5) {
                            auto& p = particles.emplace_back(C, 1.0, 10.0, globals::opt.geo_eps);
                            p.current_state.setTranslation(x * Eigen::Vector3d({x_spacing, 0, 0}) + Eigen::Vector3d({-x_offset, -y_z_offset + y_z_spacing * y + 1, -y_z_offset + y_z_spacing * y + 5}));
                            p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({globals::opt.rand_float(90)/360.0*igl::PI*2.0, globals::opt.rand_float(90)/360.0*igl::PI*2.0, globals::opt.rand_float(90)/360.0*igl::PI*2.0})));
                            p.current_state.setVelocity(Eigen::Vector3d({0, -0.1, 0}));
                        }
                    }
                }

                for (int y = 0; y < 4; y++) {
                    for (int z = 0; z < 40; z++) {
                        for (int x = 0; x < scenario_size; x++) {
                            if (globals::opt.rand_float(1.0) < 0.2) {
                                auto& p = particles.emplace_back(C, 1.0, 10.0, globals::opt.geo_eps);
                                p.current_state.setTranslation(x * Eigen::Vector3d({x_spacing, 0, 0}) + Eigen::Vector3d({-x_offset, 20 - y * 4+2, 28 + z * 10}));
                                p.current_state.setRotation(Delta2::common::eulerAnglesToQuaternion(Eigen::Vector3d({globals::opt.rand_float(90)/360.0*igl::PI*2.0, globals::opt.rand_float(90)/360.0*igl::PI*2.0, globals::opt.rand_float(90)/360.0*igl::PI*2.0})));
                                p.current_state.setVelocity(Eigen::Vector3d({0, -0.5, -globals::opt.rand_float(1.0)}));
                            }
                        }
                    }
                }
            }

            // const double offset = 5;
            // const double recentre = (double)scenario_size / 2.0 * offset;

            // for (int x = 0; x < scenario_size; x++) {
            //     double separation = globals::opt.rand_float(1.9) + 0.1;
            //     {
            //         auto& p = particles.emplace_back(M, 1.0, 10.0, globals::opt.geo_eps);
            //         p.current_state.setTranslation(x * Eigen::Vector3d({offset, 0, 0}) + Eigen::Vector3d({0.5 - recentre, -2, 0}));
            //         p.current_state.setVelocity({0, separation, 0});
            //     }
            //     {
            //         auto& p = particles.emplace_back(M, 1.0, 10.0, globals::opt.geo_eps);
            //         p.current_state.setTranslation(x * Eigen::Vector3d({offset, 0, 0}) + Eigen::Vector3d({-0.5 - recentre, 2, 0}));
            //         p.current_state.setVelocity({0, -separation, 0});
            //     }
            // }
        }       
    }
}

#include "stdio.h"
#include "../common/triangle.h"
#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../common/viewer.h"
// #include "../timestepping/explicit.h"
// #include "../timestepping/explicit_adaptive.h"
// #include "../timestepping/explicit_adaptive_clusters.h"
#include "../model/particle_handler.h"

#include <thread>
// #include <CL/sycl.hpp>
#include <ittnotify.h>

#include "../strategies/contact_force/contact_force_local.h"
#include "../strategies/contact_force/contact_force_iterative.h"
// #include "../strategies/contact_force/sequential_impulses_s_time.h"
#include "../strategies/contact_force/sequential_impulses.h"
#include "../strategies/friction/friction_iterative.h"
#include "../strategies/time_step_size/time_step_selection_dynamic_continuous.h"
#include "../strategies/contact_detection/contact_detection_comparison.h"
#include "../strategies/contact_detection_continuous/contact_detection_continuous_comparison.h"
#include "../strategies/PDE/PDE_explicit.h"
#include "../strategies/broad_phase/broad_phase_embree_cluster.h"

#include "../globals.h"

using namespace Delta2;

std::mutex globals::contact_draws_lock;
std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> globals::contact_draws;

void guiThread(common::AnimationViewer* view) {
    view->show();
}

int main(int argc, char *argv[]) {
    Delta2::common::Options opt;
    int opt_result = opt.fromArgs(argc, argv);
    if (opt_result > 0) {
        return opt_result;
    }
    
    std::vector<Delta2::Particle> particles;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // {
    //     Delta2::common::sphere(1.0, 9, V, F);

    //     std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    //     {
    //         auto& p = particles.emplace_back(M, 1.0, 0.4, 0.05);
    //         p.current_state.setTranslation({0.0, 0.0, 1.03});
    //         p.current_state.setVelocity({0.0, 0.0, -0.2});
    //         // p.current_state.setAngular({-10, -5, 0});
    //     }
    //     // {
    //     //     auto& p = particles.emplace_back(M, 1.0, 0.95, 0.05);
    //     //     p.current_state.setTranslation({3.0, 6.0, 1.1});
    //     //     p.current_state.setVelocity({0.0, -3.0, 0.0});
    //     // }
    // }

    {
        Delta2::common::plane(10.0, V, F);
        std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt, true));
        auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
        p.is_static = true;
    }

    Delta2::common::cube(V, F);

    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    // {
    //     auto& p = particles.emplace_back(M, 1.0, 0.4, 0.05);
    //     p.current_state.setTranslation({0.0, 0.0, 1.55});
    //     double angle = 45.0 / 180.0 * 3.14159;
    //     Eigen::Quaterniond r;
    //     r = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX())
    //         * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitY())
    //         * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
    //     p.current_state.setRotation(r);
    // }
    // {
    //     auto& p = particles.emplace_back(M, 1.0, 0.4, 0.5);
    //     p.current_state.setTranslation({0.0, 0.0, 1.53});
    //     p.current_state.setVelocity({0.0, 0.0, -0.2});
    //     // double angle = 1.0 / 180.0 * 3.14159;
    //     // Eigen::Quaterniond r;
    //     // r = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX())
    //     //     * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitY())
    //     //     * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
    //     // p.current_state.setRotation(r);

    // }
    // {
    //     auto& p = particles.emplace_back(M, 1.0, 0.4, 0.5);
    //     p.current_state.setTranslation({0.0, 0.0, 1.53});
    //     p.current_state.setVelocity({0.0, 1.0, 0.0});
    //     p.current_state.setAngular({0.0, 0.0, 1.0});
    // }
    // {
    //     auto& p = particles.emplace_back(M, 1.0, 0.95, 0.05);
    //     p.current_state.setTranslation({3.0, 0.0, 3.6});
    // }
    // {
    //     auto& p = particles.emplace_back(M, 1.0, 0.95, 0.25);
    //     p.current_state.setTranslation({0.0, 0.0, 1.5});
    //     p.current_state.setVelocity({10.0, 0.0, 0.0});
    //     // p.current_state.setAngular({0.0, 5.0, 10.0});
    //     // double angle = 135.0 / 180.0 * 3.14159;
    //     // Eigen::Quaterniond r;
    //     // r = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX())
    //     //     * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitY())
    //     //     * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
    //     // p.current_state.setRotation(r);
    // }
    // for (int y = 0; y < 1; y++) {
    //     for (int x = 5; x < 6; x++) {
    //         for (int i = 0; i < 10; i++) {
    //             auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
    //             p.current_state.setTranslation({-160 + x * 8.0, -160 + y * 40 + i * x * 0.2, 1.1 + i * 2.06});
    //         }
    //     }
    // }
    for (int i = 0; i < 10; i++) {
        auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
        p.current_state.setTranslation({0.0, i * 2 * 0.2, 1.1 + i * 2.06});
    }

    model::ParticleHandler ph(particles);

    std::thread gui;
    Delta2::common::AnimationViewer view(&particles);

    if (opt.gui) {
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> empty;
        view.recordFrame(empty);
        gui = std::thread(guiThread, &view);
    }

    strategy::ContactDetectionContinuousComparison contact_detection_continuous;
    strategy::TimeStepSelectionDynamicContinuous time_step(contact_detection_continuous, opt);
    strategy::ContactDetectionComparison contact_detection;
    strategy::FrictionIterative friction(opt);
    strategy::SequentialImpulses contact_force(friction, opt);
    strategy::PDEExplicit PDE(contact_detection, contact_force, friction, time_step, opt);
    strategy::BroadPhaseEmbreeCluster broad_phase(PDE, opt);

    // PDE.printType();

    time_step.init(ph);

    bool cont = opt.final_time > 0.0 || (opt.final_time < 0.0 && opt.num_time_steps > 0);
    int step = 0;
    while (cont) {
        printf("Step: %i\n", step);
        // std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> view_draws;
        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            globals::contact_draws.clear();
        }

        broad_phase.step(ph);

        if (opt.gui) {
            view.recordFrame(globals::contact_draws);
        }

        step++;
        if (opt.final_time > 0.0) {
            cont = false;
            for (const Delta2::Particle& p : particles) {
                double time = p.current_state.getTime();
                double final_time = opt.final_time;
                bool is_static = p.is_static;
                if (time < final_time && !is_static) {
                    cont = true;
                }
            }
        }
        else {
            cont = step < opt.num_time_steps;
        }
        
    }

    printf("================ Done ================\n");

    if (opt.gui) {
        gui.join();
    }
}

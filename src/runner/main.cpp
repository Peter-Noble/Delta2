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
Delta2::common::Options globals::opt;

void guiThread(common::AnimationViewer* view) {
    view->show();
}

int main(int argc, char *argv[]) {
    int opt_result = globals::opt.fromArgs(argc, argv);
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

    const int max_x = 8;
    const double spacing_x = 8.0;
    const int max_y = 8;
    const double spacing_y = 40.0;
    const int max_i = 10;
    const int tilt_offset = 4;

    // {
    //     Delta2::common::plane(std::max(10.0, std::max(max_x * spacing_x, max_y * spacing_y)), V, F);
    //     std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt, true));
    //     auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
    //     p.is_static = true;
    // }

    Delta2::common::cube(V, F);
    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, globals::opt));

    const double start_offset = 10;
    const double offset = 8;
    const int line = 10;
    const double vel_mult = 10;

    Eigen::Vector3d normal = {0, 1, 0};
    normal.normalize();

    for (int x = 0; x < line; x++) {
        auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
        p.current_state.setTranslation((start_offset + x * offset) * normal);
        p.current_state.setVelocity(-normal * vel_mult);
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

    // for (int y = 0; y < max_y; y++) {
    //     for (int x = 0; x < max_x; x++) {
    //         for (int i = 0; i < max_i; i++) {
    //             auto& p = particles.emplace_back(M, 1.0, 10.0, 0.25);
    //             p.current_state.setTranslation({-(max_x / 2 * spacing_x) + x * spacing_x, -(max_y / 2 * spacing_y) + y * spacing_y + i * (x + tilt_offset) * 0.2, 1.1 + i * 2.06});
    //         }
    //     }
    // }

    model::ParticleHandler ph(particles);

    std::thread gui;
    Delta2::common::AnimationViewer view(&particles);

    if (globals::opt.gui) {
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> empty;
        view.recordFrame(empty);
        gui = std::thread(guiThread, &view);
    }

    strategy::ContactDetectionContinuousComparison contact_detection_continuous;
    strategy::TimeStepSelectionDynamicContinuous time_step(contact_detection_continuous, globals::opt);
    strategy::ContactDetectionComparison contact_detection;
    strategy::FrictionIterative friction(globals::opt);
    strategy::SequentialImpulses contact_force(friction, globals::opt);
    strategy::PDEExplicit PDE(contact_detection, contact_force, friction, time_step, globals::opt);
    strategy::BroadPhaseEmbreeCluster broad_phase(PDE, globals::opt);

    // PDE.printType();

    time_step.init(ph);

    bool cont = globals::opt.final_time > 0.0 || (globals::opt.final_time < 0.0 && globals::opt.num_time_steps > 0);
    int step = 0;
    while (cont) {
        printf("Step: %i\n", step);
        // std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> view_draws;
        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            globals::contact_draws.clear();
        }

        broad_phase.step(ph);

        if (globals::opt.gui) {
            view.recordFrame(globals::contact_draws);
        }

        step++;
        if (globals::opt.final_time > 0.0) {
            cont = false;
            for (const Delta2::Particle& p : particles) {
                double time = p.current_state.getTime();
                double final_time = globals::opt.final_time;
                bool is_static = p.is_static;
                if (time < final_time && !is_static) {
                    cont = true;
                }
            }
        }
        else {
            cont = step < globals::opt.num_time_steps;
        }
        
    }

    printf("================ Done ================\n");

    if (globals::opt.gui) {
        gui.join();
    }
}

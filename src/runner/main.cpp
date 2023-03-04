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
// #include <ittnotify.h>

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

#include "../strategies/contact_force/resolve_penetrations_pbd.h"

#include "../scenarios/zero_g_bundle.h"
#include "../scenarios/towers.h"
#include "../scenarios/start_intersecting.h"
#include "../scenarios/pairs.h"

#include "../globals.h"

#include "tbb/tbb.h"

using namespace Delta2;

std::mutex globals::contact_draws_lock;
std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> globals::contact_draws;
Delta2::common::Options globals::opt;
Delta2::common::IttHandles globals::itt_handles;
Delta2::common::Logger globals::logger;

#if !defined(NOGL)
void guiThread(common::AnimationViewer* view) {
    view->show();
}  
#endif 

int main(int argc, char *argv[]) {
    globals::logger.printf(0, "Main\n");
    globals::itt_handles.disable_detailed_domain();

    if (globals::opt.threads > 0) {
        tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, globals::opt.threads);  // TODO Does this do something?
    }

    int opt_result = globals::opt.fromArgs(argc, argv);
    if (opt_result > 0) {
        return opt_result;
    }

    globals::logger.priority = globals::opt.print_priority;
    
    strategy::FrictionIterative friction(globals::opt);
    strategy::SequentialImpulses contact_force(friction, globals::opt);
    strategy::ContactDetectionContinuousComparison contact_detection_continuous;
    strategy::TimeStepSelectionDynamicContinuous time_step(contact_detection_continuous, globals::opt);
    strategy::ContactDetectionComparison contact_detection;
    strategy::PDEExplicit PDE(contact_detection, contact_force, friction, time_step, globals::opt);
    strategy::BroadPhaseEmbreeCluster broad_phase(PDE, globals::opt);

    std::vector<Delta2::Particle> particles;

    int scenario = 0;

    switch (globals::opt.scenario) {
    case 0:
        {
            scenarios::zero_g_bundle_orig(particles, contact_force);
            break;
        }
    case 1:
        {
            scenarios::towers(particles, contact_force);
            break;
        }
    case 2:
        {
            scenarios::start_intersecting(particles, contact_force);
            break;
        }
    case 3:
        {
            scenarios::pairs(particles, contact_force);
            break;
        }
    case 4:
        {
            scenarios::zero_g_bundle(particles, contact_force);
            break;
        }
    default:
        throw std::runtime_error("No scenario given");
    }

    model::ParticleHandler ph(particles);
    ph.initLast();

    #if !defined(NOGL)
    std::thread gui;
    Delta2::common::AnimationViewer view(&particles);

    if (globals::opt.gui) {
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> empty;
        view.recordFrame(empty);
        gui = std::thread(guiThread, &view);
    }
    #endif

    // PDE.printType();

    time_step.init(ph);

    bool cont = globals::opt.final_time > 0.0 || (globals::opt.final_time < 0.0 && globals::opt.num_time_steps > 0);
    int step = 0;
    while (cont) {
        globals::logger.printf(1, "Step: %i\n", step);
        // std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> view_draws;
        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            globals::contact_draws.clear();
        }

        broad_phase.step(ph);

        #if !defined(NOGL)
        if (globals::opt.gui) {
            view.recordFrame(globals::contact_draws);
        }
        #endif

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

    globals::logger.printf(0, "================ Done ================\n");

    #if !defined(NOGL)
    if (globals::opt.gui) {
        gui.join();
    }
    #endif
}

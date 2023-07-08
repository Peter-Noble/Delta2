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
#include "../strategies/contact_detection/contact_detection_hybrid.h"
#include "../strategies/contact_detection_continuous/contact_detection_continuous_comparison.h"
#include "../strategies/PDE/PDE_explicit.h"
#include "../strategies/broad_phase/broad_phase_embree_cluster.h"
#include "../common/export_d2.h"

#include "../strategies/contact_force/resolve_penetrations_pbd.h"

#include "../scenarios/scenarios.h"

#include "../globals.h"

#include "tbb/tbb.h"

using namespace Delta2;

std::mutex globals::contact_draws_lock;
std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> globals::contact_draws;
Delta2::common::Options globals::opt;
Delta2::common::IttHandles globals::itt_handles;
Delta2::common::Logger globals::logger;
Delta2::globals::Stats globals::stats;

#if !defined(NOGL)
void guiThread(common::AnimationViewer* view) {
    view->show();
}  
#endif 

int main(int argc, char *argv[]) {
    globals::logger.printf(0, "Main\n");
    globals::itt_handles.disable_detailed_domain();

    int opt_result = globals::opt.fromArgs(argc, argv);
    if (opt_result > 0) {
        return opt_result;
    }
    
    if (std::getenv("EXPORT_DIR") == nullptr) {
        throw std::runtime_error("EXPORT_DIR environment variable not set");
    }

    globals::logger.priority = globals::opt.print_priority;

    int threads = tbb::info::default_concurrency();
    if (globals::opt.threads > 0) {
        threads = globals::opt.threads;
    }
    globals::opt.threads = threads;
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, threads);  // TODO Does this do something?

    globals::logger.printf(1, "Threads: %i\n", threads);
    
    strategy::FrictionIterative friction(globals::opt);
    strategy::SequentialImpulses contact_force(friction, globals::opt);
    strategy::ContactDetectionContinuousComparison contact_detection_continuous;
    strategy::TimeStepSelectionDynamicContinuous time_step(contact_detection_continuous, globals::opt);
    strategy::ContactDetectionComparison contact_detection;
    // strategy::ContactDetectionHybrid contact_detection;
    strategy::PDEExplicit PDE(contact_detection, contact_force, friction, time_step, globals::opt);
    strategy::BroadPhaseEmbreeCluster broad_phase(PDE, globals::opt);

    std::vector<Delta2::Particle> particles;

    scenarios::make_scenario(particles, contact_force);

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

    Delta2::D2Writer export_writer(particles);

    time_step.init(ph);

    bool cont = globals::opt.final_time > 0.0 || (globals::opt.final_time < 0.0 && globals::opt.num_time_steps > 0);
    int step = 0;
    while (cont) {
        globals::logger.printf(1, "Step: %i\n", step);

        if (globals::opt.export_result) {
            Delta2::globals::logger.printf(3, "Starting capture\n");
            export_writer.capture(particles);
            Delta2::globals::logger.printf(3, "Done capture\n");
        }

        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            globals::contact_draws.clear();
        }

        for (Particle& p : particles) {
            const double energy_damp = 0.99;
            double damp_fac = std::pow(energy_damp, p.current_state.getTime() - p.last_state.getTime());
            p.current_state.setVelocity(p.current_state.getVelocity() * damp_fac);
            p.current_state.setAngular(p.current_state.getAngularMomentum() * damp_fac);
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

    if (globals::opt.export_result) {
        Delta2::globals::logger.printf(2, "Start final write\n");
        export_writer.capture(particles);
        std::string export_path = std::string(std::getenv("EXPORT_DIR"))+"/export.d2";
        Delta2::globals::logger.printf(1, "Exporting to %s\n", export_path.c_str());
        export_writer.write(particles, 30, export_path);
    }

    globals::logger.printf(0, "================ Done ================\n");

    if (globals::opt.collect_stats) {
        globals::stats.print();
    }

    #if !defined(NOGL)
    if (globals::opt.gui) {
        gui.join();
    }
    #endif
}

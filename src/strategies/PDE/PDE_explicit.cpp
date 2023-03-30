#include "PDE_explicit.h"

#include "../../collision_detection/full_tree_comparison.h"
#include "../../model/forces.h"
#include "../../strategies/contact_detection/contact_detection_comparison.h"
#include "../contact_force/resolve_penetrations_pbd.h"

#include "tbb/task_group.h"

#include "../../globals.h"

#include <mutex>

using namespace Delta2;
using namespace strategy;

PDEExplicit::PDEExplicit(ContactDetectionStrategy& contact_detection,
                         ContactForceStrategy& contact_force,
                         FrictionStrategy& friction,
                         TimeStepSelectionStrategy& time_step,
                         common::Options& opt) :
                         PDEStrategy(opt),
                         _contact_detection(contact_detection),
                         _contact_force(contact_force),
                         _friction(friction),
                         _time_step(time_step) {
}

double PDEExplicit::selectTimeStep(collision::Cluster& cluster) {
    // globals::logger.printf(3, "Selecting explicit time step\n");
    double time = _time_step.selectTimeStep(cluster);
    cluster.step_size = time;
    return time;
}

bool PDEExplicit::step(collision::Cluster& cluster, bool allow_fail) {
    // __itt_string_handle* compare_individual_pair_task = __itt_string_handle_create("Compare individual pair");
    
    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<int> counts;
    std::vector<bool> can_advance;

    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    for (Particle* p : cluster.particles) {
        forces.emplace_back(0.0, 0.0, 0.0);
        torques.emplace_back(0.0, 0.0, 0.0);
        counts.emplace_back(0);
    }

    std::vector<collision::Contact<double>> hits;

    for (Particle *p : cluster.particles)
    {
        p->last_time_step_size = cluster.step_size;
        if (!p->is_static)
        {
            p->rollBackState(cluster.min_current_time);
            p->projectFutureState(p->last_time_step_size);
        }
        // else
        // {
        //     p->current_state.setTime(min_times[cluster_i]);
        //     p->future_state = p->current_state;
        //     p->future_state.setTime(p->last_time_step_size);
        // }
    }

    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    std::mutex lock;

    tbb::task_group task_group;

    bool failed_before_force_solve = false;

    int failed_a_id, failed_b_id;

    // #pragma omp parallel
    // {
    //     #pragma omp single
    //     {
    //         #pragma omp taskloop
    for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
    {
        task_group.run([&, b_i] {
            if (!failed_before_force_solve) {
                // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, compare_individual_pair_task);
                collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                int a_id = cluster.particles.getLocalID(b.A.first); // geo id
                int b_id = cluster.particles.getLocalID(b.B.first);

                Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], _contact_detection);
                std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

                for (collision::Contact<double> &c : Cs)
                {
                    if ((c.A - c.B).norm() < 1e-8) {
                        failed_before_force_solve = true;
                        failed_a_id = a_id;
                        failed_b_id = b_id;
                    }
                }
                
                std::lock_guard<std::mutex> guard(lock);
                for (collision::Contact<double> &c : filtered)
                {
                    hits.push_back(c);
                }

                // __itt_task_end(globals::itt_handles.detailed_domain);
            }
        });
    }
    //     }
    // }

    task_group.wait();
    
    bool is_step_intersection_resolved = false;

    if (failed_before_force_solve) {
        if (cluster.step_size < 1e-8) {
            globals::logger.printf(3, "%i: Failed before contact solve and resolving\n", cluster_id);
            globals::logger.printf(3, "%i: Failed a: %i, failed b: %i\n", cluster_id, cluster.particles[failed_a_id].id, cluster.particles[failed_b_id].id);

            // {
            //     common::Viewer view;

            //     view.addParticleFuture(cluster.particles[failed_a_id]);
            //     view.addParticleFuture(cluster.particles[failed_b_id]);

            //     view.show();
            // }

            collision::resolvePenetrationsPBD(cluster, false);
            globals::logger.printf(3, "Resolved\n", cluster_id);

            is_step_intersection_resolved = true;

            // Eigen::Vector3d aT = cluster.particles[failed_a_id].future_state.getTranslation();
            // Eigen::Vector3d bT = cluster.particles[failed_b_id].future_state.getTranslation();

            // globals::logger.printf(3, "bT: %f, %f, %f\n", aT.x(), aT.y(), aT.z());
            // globals::logger.printf(3, "aT: %f, %f, %f\n", bT.x(), bT.y(), bT.z());

            // {
            //     common::Viewer view;

            //     view.addParticleFuture(cluster.particles[failed_a_id]);
            //     view.addParticleFuture(cluster.particles[failed_b_id]);

            //     view.show();
            // }

            // throw std::runtime_error(std::to_string(cluster_id) + std::string(": Failed before contact solve on tiny timestep"));
        }
        else {
            globals::logger.printf(3, "%i: Failed before contact solve and rolling back\n", cluster_id);
            return false;
        }
    }

    bool success = _contact_force.solve(cluster, hits, allow_fail);
    if (!success) {
        if (allow_fail) {
            globals::logger.printf(2, "%i: Failed after contact force solve\n", cluster_id);
            if (cluster.step_size < 1e-8) {
                // if (cluster.step_size < 1e-8) {
                //     double last_time, current_time;
                //     for (Particle* p : cluster.particles) {
                //         if (!p->is_static) {
                //             last_time = p->last_state.getTime();
                //             current_time = p->current_state.getTime();
                //         }
                //     }
                //     if (current_time - last_time < 1e-8) {
                //         globals::logger.printf(2, "%i: Resolving post solve fail\n", cluster_id);
                //         collision::resolvePenetrationsPBD(cluster, false);
                //     }
                // }

                collision::resolvePenetrationsPBD(cluster, false);
            }
            else {
                return false;
            }
        }
        else {
            collision::resolvePenetrationsPBD(cluster, true);
            // return false;
        }
    }
    
    strategy::ContactDetectionComparison contact_detection;
    tbb::task_group task_group_check_last;
    bool failed_at_last_after_step = false;
    bool failed_at_current_after_step = false;
    bool failed_at_future_after_step = false;


    int failed_id_a;
    int failed_id_b;

    if (globals::opt.check_post_solve) {
        for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
        {
            task_group_check_last.run([&, b_i] {
                // // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, compare_individual_pair_last_task);
                if (!failed_at_last_after_step && !failed_at_current_after_step && !failed_at_future_after_step) {
                    collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                    int a_id = cluster.particles.getLocalID(b.A.first); // geo id
                    int b_id = cluster.particles.getLocalID(b.B.first);

                    Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                    Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                    std::vector<collision::Contact<double>> Cs_last = collision::compareTreesFullLast<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);
                    std::vector<collision::Contact<double>> Cs_current = collision::compareTreesFullCurrent<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);
                    std::vector<collision::Contact<double>> Cs_future = collision::compareTreesFull<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);

                    for (collision::Contact<double> &c : Cs_last)
                    {
                        if ((c.A - c.B).norm() < 1e-4) {
                            failed_at_last_after_step = true;
                            failed_id_a = a_id;
                            failed_id_b = b_id;
                        }
                    }

                    for (collision::Contact<double> &c : Cs_current)
                    {
                        if ((c.A - c.B).norm() < 1e-4) {
                            failed_at_current_after_step = true;
                            failed_id_a = a_id;
                            failed_id_b = b_id;
                        }
                    }

                    for (collision::Contact<double> &c : Cs_future)
                    {
                        if ((c.A - c.B).norm() < 1e-4) {
                            failed_at_future_after_step = true;
                            failed_id_a = a_id;
                            failed_id_b = b_id;
                        }
                    }
                }
                // // __itt_task_end(globals::itt_handles.detailed_domain);
            });
        }

        task_group_check_last.wait();
        if (failed_at_last_after_step) {
            if (cluster.step_size < 1e-6) {
                globals::logger.printf(3, "%i: Failed at last after contact solve on tiny timestep", cluster_id);
                // common::Viewer view;

                // view.addParticle(cluster.particles[failed_id_a]);
                // view.addParticle(cluster.particles[failed_id_b]);

                // view.show();

                // throw std::runtime_error("Failed after contact solve on tiny timestep");
            }
            globals::logger.printf(3, "%i: Failed at last after step\n", cluster_id);
            return false;
        }

        if (failed_at_current_after_step && !is_step_intersection_resolved) {
            if (cluster.step_size < 1e-6) {
                globals::logger.printf(3, "%i: Failed at current after contact solve on tiny timestep", cluster_id);
                // common::Viewer view;

                // view.addParticle(cluster.particles[failed_id_a]);
                // view.addParticle(cluster.particles[failed_id_b]);

                // view.show();

                // throw std::runtime_error("Failed at current after contact solve on tiny timestep");
            }
            globals::logger.printf(3, "%i: Failed at current after step\n", cluster_id);
            return false;
        }

        if (failed_at_future_after_step) {
            if (cluster.step_size < 1e-6) {
                globals::logger.printf(3, "%i: Failed at future (%i, %i) after contact solve on tiny timestep\n", cluster_id, cluster.particles[failed_id_a].id, cluster.particles[failed_id_b].id);
                Eigen::Vector3d aT = cluster.particles[failed_id_a].future_state.getTranslation();
                Eigen::Vector3d bT = cluster.particles[failed_id_b].future_state.getTranslation();

                globals::logger.printf(3, "%i: (%f, %f, %f), (%f, %f, %f)\n", cluster_id, aT.x(), aT.y(), aT.z(), bT.x(), bT.y(), bT.z());
                // common::Viewer view;

                // view.addParticleFuture(cluster.particles[failed_id_a]);
                // view.addParticleFuture(cluster.particles[failed_id_b]);

                // view.show();

                // throw std::runtime_error("Failed future after contact solve on tiny timestep");
            }
            globals::logger.printf(3, "%i: Failed at future after step\n", cluster_id);
            return false;
        }
    }

    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    bool willSleep = true;
    for (Particle *p : cluster.particles)
    {
        if (!p->is_static)
        {
            willSleep &= p->last_state.isStationary() && p->current_state.isStationary() && p->future_state.isStationary();
            willSleep &= p->last_state.getTime() < p->current_state.getTime() && p->current_state.getTime() < p->future_state.getTime();

            p->last_state = p->current_state;
            p->current_state = p->future_state;
            p->projectFutureState(cluster.step_size);
            // globals::logger.printf(3, "Updating %i last: %.4f, current: %.4f, future: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());

            assert(p->last_state.getTime() < p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
            // globals::logger.printf(3, "Post integration state of %i: last time %f, time %f, future time %f, sleep from %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), p->sleep_candidate_time);
        }
    }

    for (Particle* p : cluster.particles)
    {
        p->setSleeping(willSleep);
    }

    return true;
}

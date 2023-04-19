#include "broad_phase_embree_cluster.h"

#include "../../model/particle.h"
#include "../../collision_detection/broad_phase_embree.h"
#include "../../collision_detection/separate_clusters.h"
#include "../../collision_detection/full_tree_comparison.h"
#include "../../common/viewer.h"
#include "../../strategies/contact_detection/contact_detection_comparison.h"
#include "../../strategies/contact_force/resolve_penetrations_pbd.h"

// #include <ittnotify.h>
// #include "tbb/task_group.h"
#include "tbb/tbb.h"

#include "../../globals.h"

using namespace Delta2;
using namespace strategy;

BroadPhaseEmbreeCluster::BroadPhaseEmbreeCluster(PDEStrategy& local_pde, common::Options& opt) :
                                                 BroadPhaseStrategy(local_pde, opt) {
    
}

void BroadPhaseEmbreeCluster::stepRecursive(Delta2::collision::Cluster& cluster, bool first_call, int depth) {
    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    if (cluster.step_size < 1e-5) {
        globals::logger.printf(3, "Small step size\n");
    }

    std::string step_size_string = std::to_string(cluster_id) + ": Step size: " + std::to_string(cluster.step_size) + ", " + (first_call ? "First" : "Following") + ". Cluster of " + std::to_string(cluster.particles.size());
    step_size_string += ". {";
    for (Particle* p : cluster.particles) {
        step_size_string += std::to_string(p->id) + ", ";
    }
    step_size_string += "}\n";

    globals::logger.printf(2, step_size_string.c_str());
    
    // assert(cluster.step_size > 0.0);
    // __itt_string_handle* step_recursive_task = __itt_string_handle_create("Step recursive task");
    // __itt_string_handle* local_pde_task = __itt_string_handle_create("Local pde step");
    // __itt_string_handle* compare_individual_pair_current_task = __itt_string_handle_create("Compare individual pair current");
    // __itt_string_handle* compare_individual_pair_last_task = __itt_string_handle_create("Compare individual pair last");
    // __itt_string_handle* roll_back_to_valid_current_state_task = __itt_string_handle_create("Roll back to valid current state");

    // __itt_string_handle* step_rec_after_first_valid_task = __itt_string_handle_create("Step rec after first valid");
    // __itt_string_handle* step_rec_after_first_valid_sep_task = __itt_string_handle_create("Step rec after first valid separated");
    // __itt_string_handle* step_rec_after_first_invalid_task = __itt_string_handle_create("Step rec after first invalid");
    // __itt_string_handle* step_rec_div_task = __itt_string_handle_create("Step rec div");
    // __itt_string_handle* step_rec_div_sep_task = __itt_string_handle_create("Step rec div separated");
    
    double orig_step_size = cluster.step_size;
    double time = std::numeric_limits<double>::infinity();
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            if (std::isinf(time)) {
                time = p->current_state.getTime();
            }
            else {
                assert(std::fabs(p->current_state.getTime() - time) < 1e-5);
            }
        }
    }

    double target_time = time + cluster.step_size;

    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_recursive_task);

    double start_step = cluster.step_size;
    bool success = false;
    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         globals::logger.printf(3, "stepRecursive  pre: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }

    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, local_pde_task);
    const int max_recursive_step_depth = 6;
    bool forced_advance = depth >= max_recursive_step_depth && cluster.step_size < 1e-8;
    globals::logger.printf(3, "%i: Step recursive depth: %i\n", cluster_id, depth);
    if (forced_advance) {
        globals::logger.printf(3, "%i: Forced advance\n", cluster_id);
        collision::resolvePenetrationsPBD(cluster, true);
    }
    success = _local_pde.step(cluster, forced_advance);
    if (forced_advance) {
        collision::resolvePenetrationsPBD(cluster, true);
    }
    // __itt_task_end(globals::itt_handles.detailed_domain);

    strategy::ContactDetectionComparison contact_detection;

    if (success) {
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                assert(std::fabs(p->current_state.getTime() - (time + cluster.step_size)) < 1e-5);
            }
        }
        globals::logger.printf(2, "local pde claims correct advance\n");
    }
    else {
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                assert(std::fabs(p->current_state.getTime() - time) < 1e-5);
            }
        }
    }

    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         globals::logger.printf(3, "stepRecursive post: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }
    if (!success && !forced_advance) {
        if (!forced_advance) {
            globals::logger.printf(2, "%i: Failed but didn't force advance\n", cluster_id);
        }
        globals::logger.printf(3, "%i: Fail\n", cluster_id);
        // Interpolate current_new back to 50% between last_start and current_start
        // Set future_new to current_start
        // Step forward from there to the current_start (recursive call)
        // Store current_new (to use in last_final)
        // Set timestep size to step_size_start / 2.0 and project new future state
        // Step forward (recursive call)
        // Step forward (recursive call)
        // Set last_final to the stored current_new and each of the particles last_time_step_size to min(last_time_step_size, step_size_start / 2.0)
        
        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, roll_back_to_valid_current_state_task);

        bool failed_at_current = true;
        bool failed_at_last = false;
        bool first_rollback = true;

        double second_step_size = 0.0;
        

        tbb::task_group task_group_check_last;

        for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
        {
            task_group_check_last.run([&, b_i] {
                // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, compare_individual_pair_last_task);
                if (!failed_at_last) {
                    collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                    int a_id = cluster.particles.getLocalID(b.A.first); // geo id
                    int b_id = cluster.particles.getLocalID(b.B.first);

                    Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                    Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                    std::vector<collision::Contact<double>> Cs = collision::compareTreesFullLast<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);

                    for (collision::Contact<double> &c : Cs)
                    {
                        if ((c.A - c.B).norm() < 1e-4) {
                            failed_at_last = true;
                        }
                    }
                }
                // __itt_task_end(globals::itt_handles.detailed_domain);
            });
        }

        task_group_check_last.wait();

        std::vector<double> tmp_orig_times;
        for (Particle* p : cluster.particles) {
            tmp_orig_times.push_back(p->current_state.getTime());
        }

        bool rolling_back = false;
        bool rolling_back_to_start = false;

        int c = 0;
        const int max_c = 4;
        while (failed_at_current) {
            c++;
            if (c > max_c) {
                globals::logger.printf(3, "%i: Can't find valid last state\n", cluster_id);
                // #ifndef NDEBUG
                // throw std::runtime_error("Can't find valid last state");
                // #endif
                break;
            }
            failed_at_current = false;

            double new_min_time = 0.0;
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    // assert(p->current_state.getTime() - p->last_state.getTime() > 1e-6);
                    new_min_time = std::max(new_min_time, common::lerp(p->last_state.getTime(), p->current_state.getTime(), std::pow(0.5, c)));
                    if (new_min_time - p->last_state.getTime() < 1e-5) {
                        c = max_c;
                    }
                }
            }

            cluster.min_current_time = new_min_time; // in the next stepRecursive call all the particles will get rolled back to this time
            // globals::logger.printf(3, "Rolling back to %f\n", cluster.min_current_time);

            double step_size;
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    p->future_state = p->current_state;
                    if (c == max_c) {
                        p->current_state = p->last_state;
                        rolling_back = true;
                        rolling_back_to_start = true;
                        cluster.min_current_time = p->current_state.getTime();
                    }
                    else {
                        p->current_state = p->current_state.interpolate(p->last_state, new_min_time);
                        rolling_back = true;
                    }
                    // cluster.min_current_time = std::min(cluster.min_current_time, p->current_state.getTime());
                    step_size = p->future_state.getTime() - p->current_state.getTime();
                    // assert(step_size > 0.0);
                }
            }

            tbb::task_group task_group;

            common::Edge failed_at(Eigen::Vector3d({0,0,0}), Eigen::Vector3d({0,0,0}));

            for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
            {
                task_group.run([&, b_i] {
                    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, compare_individual_pair_current_task);
                    if (!failed_at_current) {
                        collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                        int a_id = cluster.particles.getLocalID(b.A.first); // geo id
                        int b_id = cluster.particles.getLocalID(b.B.first);

                        Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                        Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                        std::vector<collision::Contact<double>> Cs = collision::compareTreesFullCurrent<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);
                        // std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

                        // std::vector<collision::Contact<double>> Csl = collision::compareTreesFullLast<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);

                        // bool failed_at_last = false;
                        // for (collision::Contact<double> &c : Csl)
                        // {
                        //     if ((c.A - c.B).norm() < 1e-4) {
                        //         failed_at_current = failed_at_last;
                        //     }
                        // }

                        for (collision::Contact<double> &c : Cs)
                        {
                            // if ((c.A - c.B).norm() < 1e-4 && failed_at_last) {
                            if ((c.A - c.B).norm() < 1e-4) {
                                failed_at_current = true;
                                failed_at = common::Edge(c.A, c.B);
                            }
                        }
                    }
                    // __itt_task_end(globals::itt_handles.detailed_domain);
                });
            }

            task_group.wait();
            cluster.step_size = step_size;

            // if (failed_at_current) {
            //     common::Viewer view;
            //     for (Particle* p : cluster.particles) {
            //         view.addParticle(*p);
            //     }
            //     view.addEdge(failed_at);
            //     view.show();
            // }

            if (!first_rollback) {
            }

            first_rollback = false;

            // if (!failed_at_current) {
            //     globals::logger.printf(3, "Rolling back to %.4f.  First time step: %.4f, second time step: %.4f.  Original current: %.4f\n", new_min_time, step_size, second_step_size);
            // }
            if (failed_at_current) {
                second_step_size += step_size;
            }
        }
        // __itt_task_end(globals::itt_handles.detailed_domain);

        if (rolling_back && !rolling_back_to_start) {
            globals::logger.printf(3, "%i: Rolling back at %f. Cluster of %i\n", cluster_id, cluster.min_current_time, cluster.particles.size());
        }
        if (rolling_back_to_start) {
            globals::logger.printf(3, "%i: Rolling back to last at %f. Cluster of %i\n", cluster_id, cluster.min_current_time, cluster.particles.size());
        }

        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            if (!cluster.particles[p_i].is_static) {
                double current_time = cluster.particles[p_i].current_state.getTime();

                second_step_size = time - current_time - cluster.step_size;

                double advanced = cluster.step_size + second_step_size;
                double difference = current_time + advanced - time;
                assert(difference < 1e-5);
            }
        }
        for (Particle* p : cluster.particles) {
            assert(p->current_state.getTime() + (cluster.step_size + second_step_size) - time < 1e-5);
        }

        // second_step_size -= cluster.step_size;

        // globals::logger.printf(3, "Before after first valid\n");
        // for (Particle* p : cluster.particles) {
        //     if (!p->is_static) {
        //         globals::logger.printf(3, "Particle %i last: %.4f, current: %.4f, future: %.4f, cluster step: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), cluster.step_size);
        //     }
        // }

        // std::vector<collision::Cluster> after_first_valid_separated_clusters = separateClusterByTimestep(cluster);
        std::vector<collision::Cluster> after_first_valid_separated_clusters = {};

        if (after_first_valid_separated_clusters.size() > 1) {
            tbb::task_group after_valid_task_group;

            for (int sub_ci = 0; sub_ci < after_first_valid_separated_clusters.size(); sub_ci++) {
                after_valid_task_group.run([&, sub_ci] {
                    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_rec_after_first_valid_sep_task);
                    stepRecursive(after_first_valid_separated_clusters[sub_ci], false, depth + 1);
                    // __itt_task_end(globals::itt_handles.detailed_domain);
                });
            }

            after_valid_task_group.wait();
        }
        else {
            // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_rec_after_first_valid_task);
            stepRecursive(cluster, false, depth + 1);
            // __itt_task_end(globals::itt_handles.detailed_domain);
        }


        cluster.step_size = second_step_size;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                cluster.min_current_time = p->current_state.getTime();
                p->projectFutureState(cluster.step_size);

                assert(std::fabs(p->current_state.getTime() + cluster.step_size - time) < 1e-5);
            }
        }

        if (second_step_size > 0.0) {
            // globals::logger.printf(3, "Before after first invalid\n");
            // for (Particle* p : cluster.particles) {
            //     if (!p->is_static) {
            //         globals::logger.printf(3, "Particle %i last: %.4f, current: %.4f, future: %.4f, cluster step: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), cluster.step_size);
            //     }
            // }

            // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_rec_after_first_invalid_task);
            stepRecursive(cluster, false, depth + 1);
            // __itt_task_end(globals::itt_handles.detailed_domain);
            
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    // Should be back to the same current time
                    assert(std::fabs(p->current_state.getTime() - time) < 1e-5);
                }
            }
        }

        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                // Should be back to the same current time
                assert(std::fabs(p->current_state.getTime() - time) < 1e-5);
            }
        }

        std::vector<State> last_state_final;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                last_state_final.push_back(p->current_state);
            }
        }

        int sub_steps = 2;

        cluster.step_size = start_step / (double) sub_steps;

        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                cluster.min_current_time = p->current_state.getTime();
                p->projectFutureState(cluster.step_size);

                assert(std::fabs((p->current_state.getTime() + sub_steps * cluster.step_size) - (time + orig_step_size)) < 1e-5);
            }
        }

        // globals::logger.printf(3, "Before sub stepping\n");
        // for (Particle* p : cluster.particles) {
        //     if (!p->is_static) {
        //         globals::logger.printf(3, "Particle %i last: %.4f, current: %.4f, future: %.4f, cluster step: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), cluster.step_size);
        //     }
        // }

        // double time = std::numeric_limits<double>::infinity();

        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_rec_div_task);
        for (int s = 0; s < sub_steps; s++) {
            // std::vector<collision::Cluster> sub_step_separated_clusters = separateClusterByTimestep(cluster, s * cluster.step_size);
            std::vector<collision::Cluster> sub_step_separated_clusters = {};
            if (sub_step_separated_clusters.size() > 1) {
                tbb::task_group sub_step_task_group;
                globals::logger.printf(3, "%i: Separating sub step clusters into %i\n", cluster_id, sub_step_separated_clusters.size());

                for (int sub_ci = 0; sub_ci < sub_step_separated_clusters.size(); sub_ci++) {
                    sub_step_task_group.run([&, sub_ci] {

                    for (Particle* p : sub_step_separated_clusters[sub_ci].particles) {
                        if (!p->is_static) {
                            assert(std::fabs((p->current_state.getTime() + (sub_steps - s) * cluster.step_size) - (time + orig_step_size)) < 1e-5);
                        }
                    }

                    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, step_rec_div_sep_task);
                    stepRecursive(sub_step_separated_clusters[sub_ci], false, depth + 1);
                    // __itt_task_end(globals::itt_handles.detailed_domain);
                
                    for (Particle* p : sub_step_separated_clusters[sub_ci].particles) {
                        if (!p->is_static) {
                            assert(std::fabs((p->current_state.getTime() + (sub_steps - s - 1) * cluster.step_size) - (time + orig_step_size)) < 1e-5);
                        }
                    }

                    });
                }

                sub_step_task_group.wait();
            }
            else {
                double target_time_tmp;
                double t0_tmp;
                for (Particle* p : cluster.particles) {
                    if (!p->is_static) {
                        assert(std::fabs((p->current_state.getTime() + (sub_steps - s) * cluster.step_size) - (time + orig_step_size)) < 1e-5);
                        t0_tmp = p->current_state.getTime();
                        target_time_tmp = t0_tmp + cluster.step_size;
                    }
                }

                stepRecursive(cluster, false, depth + 1);

                for (Particle* p : cluster.particles) {
                    if (!p->is_static) {
                        assert(std::fabs(p->current_state.getTime() - (t0_tmp + cluster.step_size) < 1e-5));
                        assert(std::fabs(p->current_state.getTime() - target_time_tmp) < 1e-5);
                        assert(std::fabs((p->current_state.getTime() + (sub_steps - s - 1) * cluster.step_size) - (time + orig_step_size)) < 1e-5);
                    }
                }
            }

            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    cluster.min_current_time = p->current_state.getTime();
                    // globals::logger.printf(3, "Step 3 %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());

                    assert(std::fabs((p->current_state.getTime() + (sub_steps - s - 1) * cluster.step_size) - (time + orig_step_size)) < 1e-5);
                }
            }
        }
        // __itt_task_end(globals::itt_handles.detailed_domain);

        // for (Particle* p : cluster.particles) {
        //     if (!p->is_static) {
        //         if (!(std::fabs(p->current_state.getTime() - (time + orig_step_size)) < 1e-5)) {
        //             globals::logger.printf(3, "Particle: %i, current: %.4f, target: %.4f\n", p->id, p->current_state.getTime(), time + orig_step_size);
        //         }
        //         assert(std::fabs(p->current_state.getTime() - (time + orig_step_size)) < 1e-5);
        //     }
        // }

        int i = 0;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                p->last_state = last_state_final[i++];
                p->last_time_step_size = std::min(p->last_time_step_size, start_step / 2.0);
            }
        }        
    } 
    // else {
    //     globals::logger.printf(3, "Success\n");
    // }

    // cluster.step_size = start_step;

    // globals::logger.printf(3, "Cluster step size: %f\n", cluster.step_size);
    double final_time = std::numeric_limits<double>::infinity();
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            assert(std::fabs(p->current_state.getTime() - target_time) < 1e-5);

            if (std::isinf(final_time)) {
                final_time = p->current_state.getTime();
            }
            else {
                assert(std::fabs(p->current_state.getTime() - final_time) < 1e-5);
            }
        }
    }

    cluster.step_size = orig_step_size;

    // __itt_task_end(globals::itt_handles.detailed_domain);
}

void BroadPhaseEmbreeCluster::step(model::ParticleHandler& particles) {
    __itt_string_handle* broad_phase_step_task = __itt_string_handle_create("Broad phase step");
    __itt_string_handle* select_time_step_task = __itt_string_handle_create("Select time phase");
    __itt_string_handle* separate_collision_clusters_task = __itt_string_handle_create("Separate collision cluster phase");
    __itt_string_handle* select_and_sort_task = __itt_string_handle_create("Select and sort clusters phase");
    __itt_string_handle* outer_step_task = __itt_string_handle_create("Outer step phase");
    __itt_string_handle* advance_cluster_task = __itt_string_handle_create("Advance cluster");
    __itt_string_handle* embree_task = __itt_string_handle_create("Embree");

    __itt_task_begin(globals::itt_handles.step_domain, __itt_null, __itt_null, broad_phase_step_task);

    for (Particle* p : particles)
    {
        p->last_time_step_size *= 2;
        p->last_time_step_size = std::min(p->last_time_step_size, (double) _opt.time_step_size);
        p->projectFutureState(p->last_time_step_size);
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, embree_task);
    collision::BroadPhaseCollisions B;
    B = collision::broadPhaseEmbree(particles);
    __itt_task_end(globals::itt_handles.phases_domain);


    __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, separate_collision_clusters_task);
    std::vector<collision::Cluster> clusters = collision::separateCollisionClusters(B, particles);
    // globals::logger.printf(3, "%i clusters found from separateCollisionClusters\n", clusters.size());
    __itt_task_end(globals::itt_handles.phases_domain);
    
    double min_final_time = std::numeric_limits<double>::infinity();

    for (collision::Cluster& cluster : clusters) {
        if (!cluster.is_static) {
            double cluster_min_time = cluster.min_current_time + cluster.step_size;
            min_final_time = std::min(min_final_time, cluster_min_time);
        }
    }

    // globals::logger.printf(3, "Clusters: %i\n", clusters.size());

    std::mutex fine_min_final_time_lock;
    double fine_min_final_time = std::numeric_limits<double>::infinity();

    // globals::logger.printf(3, "Select time step sizes\n");
    tbb::task_group task_group;

    std::vector<collision::Cluster> clusters_tmp;
    clusters_tmp.reserve(clusters.size());

    __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, select_time_step_task);

    int total_clusters = clusters.size();
    for (int cluster_i = 0; cluster_i < total_clusters; cluster_i++)
    {
        // globals::logger.printf(3, "Considering cluster %i\n", cluster_i);
        collision::Cluster& cl = clusters[cluster_i];
        // clusters.pop_back();

        auto it = std::find_if(last_step_clusters.begin(), last_step_clusters.end(), [&cl](const collision::Cluster& last) {return !last.hasAdvanced(cl);});        
        if (it != last_step_clusters.end()) {
            // This cluster was the same as it was last iteration (and it didn't advance) so just keep the same step size recommendation.
            if (cl.step_size >= 0) {
                // If step_size has been set to -1 then this cluster has been seen but not advanced before
                cl.step_size = it->step_size;
                fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + it->step_size);
                // globals::logger.printf(1, "%i skipping time step selection because time step has already been picked for this cluster since last modification\n", cluster_i);
                // clusters_tmp.push_back(cl);
                continue;
            }
        }

        if (cl.min_current_time > min_final_time)
        {
            // globals::logger.printf(1, "%i skipping time step selection because min current time ahead of min final time\n", cluster_i);
            // clusters_tmp.push_back(cl);
            // This cluster can't advance because there is another cluster that is too far behind
            continue;
        }

        if (cl.is_static)
        {
            // globals::logger.printf(1, "%i skipping time step selection because static cluster\n", cluster_i);
            // This cluster shouldn't be advanced as it's asleep or static
            // clusters_tmp.push_back(cl);
            continue;
        }

        if (cl.sleeping) {
            // globals::logger.printf(1, "%i skipping time step selection because sleeping cluster\n", cluster_i);
            fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + _opt.time_step_size);
            cl.step_size = _opt.time_step_size;
            // clusters_tmp.push_back(cl);
            continue;
        }
        
        task_group.run([&, cluster_i] {  // TODO this causes segfaults when used with cluster separation for some reason
            collision::Cluster& cluster = clusters[cluster_i];

            double step = _local_pde.selectTimeStep(cluster);

            std::string step_size_string = "Selected time step size: " + std::to_string(step) + " for cluster of " + std::to_string(cluster.particles.size());
            step_size_string += ". {";
            for (Particle* p : cluster.particles) {
                step_size_string += std::to_string(p->id) + ", ";
            }
            step_size_string += "}\n";

            globals::logger.printf(2, step_size_string.c_str());
            // std::vector<collision::Cluster> sub_clusters = collision::separateClusterByTimestep(clusters[cluster_i]);
            // std::vector<collision::Cluster> sub_clusters = {clusters[cluster_i]};

            {
                std::lock_guard guard(fine_min_final_time_lock);
                fine_min_final_time = std::min(fine_min_final_time, clusters[cluster_i].min_current_time + step);

                // globals::logger.printf(1, "fine_min_final_time: %f\n", fine_min_final_time);
                // globals::logger.printf(1, "Selected time: %f for cluster %i\n", step, cluster_i);

                // if (sub_clusters.size() > 1) {
                //     int up_to_cluster = clusters.size();
                //     for (int sub_i = 0; sub_i < sub_clusters.size(); sub_i++) {
                //         collision::Cluster sub_cl = sub_clusters[sub_i];
                //         for (Particle* p : sub_cl.particles) {
                //             p->cluster_id = up_to_cluster + sub_i;
                //         }
                //         // if (sub_i == 0) {
                //         //     clusters[cluster_i] = sub_cl;
                //         // }
                //         // else {
                //         //     clusters.push_back(sub_cl);
                //         // }
                //         // clusters_tmp.push_back(sub_cl);
                //     }

                //     globals::logger.printf(3, "Splitting into %i\n", sub_clusters.size());
                // }
            }
        });
    }
    task_group.wait();
    __itt_task_end(globals::itt_handles.phases_domain);

    if (!globals::opt.local_ts) {
        double accum_step = globals::opt.time_step_size;

        for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++) {
            accum_step = std::min(accum_step, clusters[cluster_i].step_size);
        }

        for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++) {
            clusters[cluster_i].step_size = accum_step;
            for (Particle* p : clusters[cluster_i].particles) {
                p->projectFutureState(accum_step);
                // globals::logger.printf(1, "Future time for particle %i: %f\n", p->id, p->future_state.getTime());
            }
        }
    }

    // clusters.clear();
    // clusters.swap(clusters_tmp);
    clusters_tmp.clear();

    for (Particle* p : particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    // Splitting into two loops makes this point act like a sync point.
    //    Advantage - Only clusters that can advance are advanced.
    //    Disadvantage - We have to wait for all time step selections to be done.

    __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, select_and_sort_task);

    std::vector<int> can_advance;

    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        if (globals::opt.local_ts) {
            if (clusters[cluster_i].min_current_time > fine_min_final_time)
            {
                continue;
            }

            if (clusters[cluster_i].is_static)
            {
                continue;
            }
        }

        if (clusters[cluster_i].sleeping) {
            _local_pde.stepSleeping(clusters[cluster_i]);
        }
        else {
            can_advance.push_back(cluster_i);
        }
    }

    if (globals::opt.local_ts) {
        std::sort(can_advance.begin(), can_advance.end(), [&clusters](int a, int b) {
            return clusters[a].min_current_time < clusters[b].min_current_time;
        });
    }

    __itt_task_end(globals::itt_handles.phases_domain);

    tbb::task_group step_task_group;

    int cluster_limit = can_advance.size();
    if (globals::opt.limit_clusters_per_step && globals::opt.local_ts) {
        int threads = globals::opt.threads > 0 ? globals::opt.threads : tbb::info::default_concurrency();
        
        cluster_limit = std::min((int)can_advance.size(), std::max(threads * 2, (int)can_advance.size() - threads));
    }

    {
        double min_current = std::numeric_limits<double>::infinity();
        for (Particle* p : particles) {
            if (!p->is_static) {
                assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
                min_current = std::min(min_current, p->current_state.getTime());
                // globals::logger.printf(3, "Particle %i last: %f, current: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
            }
        }
        for (Particle* p : particles) {
            if (!p->is_static) {
                assert(p->last_state.getTime() <= min_current);
            }
        }
    }

    __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, outer_step_task);
    // globals::logger.printf(3, "Advance\n");
    // for (int ca = 0; ca < std::min((int)can_advance.size(), num_tasks); ca++)
    for (int ca = 0; ca < cluster_limit; ca++)
    {
        int cluster_i = can_advance[ca];
        step_task_group.run([&, cluster_i] {
            __itt_task_begin(globals::itt_handles.phases_domain, __itt_null, __itt_null, advance_cluster_task);

            // globals::logger.printf(1, "cluster_i: %i, min_current_time: %f, step_size: %f\n", cluster_i, clusters[cluster_i].min_current_time, clusters[cluster_i].step_size);
            assert(clusters[cluster_i].min_current_time + clusters[cluster_i].step_size >= fine_min_final_time);
            for (Particle* p : clusters[cluster_i].particles) {
                if (!p->is_static) {
                    // globals::logger.printf(1, "Pre  step particle %i last: %.4f, current: %.4f, future: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
                }
            }

            stepRecursive(clusters[cluster_i], true);
            for (Particle* p : clusters[cluster_i].particles) {
                if (!p->is_static) {
                    // globals::logger.printf(1, "Post step particle %i last: %.4f, current: %.4f, future: %.4f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
                }
            }

            for (Particle* p : clusters[cluster_i].particles) {
                if (!p->is_static) {
                    if (p->current_state.getTime() < fine_min_final_time) {
                        if (p->current_state.getTime() + 1e-8 > fine_min_final_time) {
                            p->current_state.setTime(fine_min_final_time);
                        }
                        else {
                            assert(false);
                        }
                    }
                }
            }
            __itt_task_end(globals::itt_handles.phases_domain);
        });
    }
    step_task_group.wait();

    __itt_task_end(globals::itt_handles.phases_domain);

    if (!globals::opt.local_ts) {
        for (Particle* p : particles) {
            if (!p->is_static && std::fabs(fine_min_final_time - p->current_state.getTime()) > 1e-6) {
                throw std::runtime_error("Time is different but no local ts");
            }
        }
    }

    for (int cluster_i = total_clusters - 1; cluster_i >= 0; cluster_i--)
    {
        collision::Cluster& cl = clusters[cluster_i];

        if (cl.min_current_time > min_final_time || cl.is_static || cl.sleeping) {
            // Remove this cluster in the next cache (stops clusters that haven't had a time step selected from being used)
            clusters.erase(clusters.begin() + cluster_i);
        }
    }

    last_step_clusters.swap(clusters);

    __itt_task_end(globals::itt_handles.step_domain);
}

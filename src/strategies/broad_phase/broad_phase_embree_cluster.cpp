#include "broad_phase_embree_cluster.h"

#include "../../model/particle.h"
#include "../../collision_detection/broad_phase_embree.h"
#include "../../collision_detection/separate_clusters.h"
#include "../../collision_detection/full_tree_comparison.h"
#include "../../common/viewer.h"
#include "../../strategies/contact_detection/contact_detection_comparison.h"

#include <ittnotify.h>
#include "tbb/task_group.h"

using namespace Delta2;
using namespace strategy;

BroadPhaseEmbreeCluster::BroadPhaseEmbreeCluster(PDEStrategy& local_pde, common::Options& opt) :
                                                 BroadPhaseStrategy(local_pde, opt) {
    
}

void BroadPhaseEmbreeCluster::stepRecursive(Delta2::collision::Cluster& cluster, bool first_call, int depth) {
    printf("Step size: %f, %s\n", cluster.step_size, first_call ? "First" : "Following");
    assert(cluster.step_size > 0.0);
    __itt_domain* domain = __itt_domain_create("My Domain");
    __itt_string_handle* step_recursive_task = __itt_string_handle_create("Step recursive task");
    __itt_string_handle* local_pde_task = __itt_string_handle_create("Local pde step");
    __itt_string_handle* compare_individual_pair_current_task = __itt_string_handle_create("Compare individual pair current");
    __itt_string_handle* roll_back_to_valid_current_state_task = __itt_string_handle_create("Roll back to valid current state");

    __itt_string_handle* step_rec_after_first_valid_task = __itt_string_handle_create("Step rec after first valid");
    __itt_string_handle* step_rec_after_first_invalid_task = __itt_string_handle_create("Step rec after first invalid");
    __itt_string_handle* step_rec_div_task = __itt_string_handle_create("Step rec div");
    
    __itt_task_begin(domain, __itt_null, __itt_null, step_recursive_task);

    double start_step = cluster.step_size;
    bool success = false;
    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         printf("stepRecursive  pre: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }

    __itt_task_begin(domain, __itt_null, __itt_null, local_pde_task);
    success = _local_pde.step(cluster, depth < 2);
    __itt_task_end(domain);

    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         printf("stepRecursive post: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }
    if (!success) {
        // printf("Fail\n");
        // Interpolate current_new back to 50% between last_start and current_start
        // Set future_new to current_start
        // Step forward from there to the current_start (recursive call)
        // Store current_new (to use in last_final)
        // Set timestep size to step_size_start / 2.0 and project new future state
        // Step forward (recursive call)
        // Step forward (recursive call)
        // Set last_final to the stored current_new and each of the particles last_time_step_size to min(last_time_step_size, step_size_start / 2.0)
        
        __itt_task_begin(domain, __itt_null, __itt_null, roll_back_to_valid_current_state_task);
        strategy::ContactDetectionComparison contact_detection;

        bool failed_at_current = true;
        bool first_rollback = true;

        double second_step_size = 0.0;
        
        int c = 0;
        while (failed_at_current) {
            c++;
            if (c > 10) {
                printf("Can't find valid last state\n");
                break;
            }
            failed_at_current = false;

            double new_min_time = 0.0;
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    new_min_time = std::max(new_min_time, (p->last_state.getTime() + p->current_state.getTime()) / 2.0);
                }
            }

            cluster.min_current_time = new_min_time; // in the next stepRecursive call all the particles will get rolled back to this time
            // printf("Rolling back to %f\n", cluster.min_current_time);

            double step_size;
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    p->future_state = p->current_state;
                    p->current_state = p->current_state.interpolate(p->last_state, new_min_time);
                    // cluster.min_current_time = std::min(cluster.min_current_time, p->current_state.getTime());
                    step_size = p->future_state.getTime() - new_min_time;
                    assert(step_size > 0.0);
                }
            }

            tbb::task_group task_group;

            common::Edge failed_at(Eigen::Vector3d({0,0,0}), Eigen::Vector3d({0,0,0}));

            for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
            {
                task_group.run([&, b_i] {
                    __itt_task_begin(domain, __itt_null, __itt_null, compare_individual_pair_current_task);
                    if (!failed_at_current) {
                        collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                        int a_id = cluster.particles.getLocalID(b.first.first); // geo id
                        int b_id = cluster.particles.getLocalID(b.second.first);

                        Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                        Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                        std::vector<collision::Contact<double>> Cs = collision::compareTreesFullCurrent<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);
                        // std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

                        for (collision::Contact<double> &c : Cs)
                        {
                            if ((c.A - c.B).norm() < 1e-4) {
                                failed_at_current = true;
                                failed_at = common::Edge(c.A, c.B);
                            }
                        }
                    }
                    __itt_task_end(domain);
                });
            }

            task_group.wait();
            cluster.step_size = step_size;

            if (failed_at_current) {
                common::Viewer view;
                for (Particle* p : cluster.particles) {
                    view.addParticle(*p);
                }
                view.addEdge(failed_at);
                view.show();
            }

            if (!first_rollback) {
                second_step_size += step_size;
            }

            first_rollback = false;
        }
        __itt_task_end(domain);

        // second_step_size -= cluster.step_size;

        __itt_task_begin(domain, __itt_null, __itt_null, step_rec_after_first_valid_task);
        stepRecursive(cluster, false, depth + 1);
        __itt_task_end(domain);

        cluster.step_size = second_step_size;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                cluster.min_current_time = p->current_state.getTime();
                p->projectFutureState(cluster.step_size);
            }
        }

        if (second_step_size > 0.0) {
            __itt_task_begin(domain, __itt_null, __itt_null, step_rec_after_first_invalid_task);
            stepRecursive(cluster, false, depth + 1);
            __itt_task_end(domain);
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
            }
        }

        __itt_task_begin(domain, __itt_null, __itt_null, step_rec_div_task);
        for (int s = 0; s < sub_steps; s++) {
            stepRecursive(cluster, false, depth + 1);
            for (Particle* p : cluster.particles) {
                if (!p->is_static) {
                    cluster.min_current_time = p->current_state.getTime();
                    // printf("Step 3 %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
                }
            }
        }
        __itt_task_end(domain);

        int i = 0;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                p->last_state = last_state_final[i++];
                p->last_time_step_size = std::min(p->last_time_step_size, start_step / 2.0);
            }
        }        
    }
    else {
        // printf("Success\n");
    }

    cluster.step_size = start_step;

    __itt_task_end(domain);
}

void BroadPhaseEmbreeCluster::step(model::ParticleHandler& particles) {
    __itt_domain* domain = __itt_domain_create("My Domain");
    __itt_string_handle* broad_phase_step_task = __itt_string_handle_create("Broad phase step");
    __itt_string_handle* select_time_step_task = __itt_string_handle_create("Select time step");
    __itt_string_handle* separate_collision_clusters_task = __itt_string_handle_create("Separate collision cluster");
    __itt_string_handle* embree_task = __itt_string_handle_create("Embree");

    __itt_task_begin(domain, __itt_null, __itt_null, broad_phase_step_task);

    for (Particle* p : particles)
    {
        p->last_time_step_size *= 2;
        p->last_time_step_size = std::min(p->last_time_step_size, (double) _opt.time_step_size);
        p->projectFutureState(p->last_time_step_size);
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    __itt_task_begin(domain, __itt_null, __itt_null, embree_task);
    collision::BroadPhaseCollisions B;
    B = collision::broadPhaseEmbree(particles);
    __itt_task_end(domain);


    __itt_task_begin(domain, __itt_null, __itt_null, separate_collision_clusters_task);
    std::vector<collision::Cluster> clusters = collision::separateCollisionClusters(B, particles);
    __itt_task_end(domain);
    
    double min_final_time = std::numeric_limits<double>::infinity();

    for (collision::Cluster& cluster : clusters) {
        if (!cluster.is_static) {
            double cluster_min_time = cluster.min_current_time + cluster.step_size;
            min_final_time = std::min(min_final_time, cluster_min_time);
        }
    }

    // printf("Clusters: %i\n", clusters.size());

    std::mutex fine_min_final_time_lock;
    double fine_min_final_time = std::numeric_limits<double>::infinity();

    // printf("Select time step sizes\n");
    tbb::task_group task_group;

    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        task_group.run([&, cluster_i] {
            collision::Cluster& cl = clusters[cluster_i];

            if (cl.min_current_time > min_final_time)
            {
                printf("%i skipping time step selection because min current time ahead of min final time\n", cluster_i);
                // This cluster can't advance because there is another cluster that is too far behind
                return;
            }

            if (cl.is_static)
            {
                // printf("%i skipping time step selection because static cluster\n", cluster_i);
                // This cluster shouldn't be advanced as it's asleep or static
                return;
            }

            if (cl.sleeping) {
                printf("%i skipping time step selection because sleeping cluster\n", cluster_i);
                fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + _opt.time_step_size);
                return;
            }

            // TODO can this test ever pass?  Does the first check for min_current_time > min_final_time catch this?
            auto it = std::find_if(last_step_clusters.begin(), last_step_clusters.end(), [&cl](const collision::Cluster& last) {return !last.hasAdvanced(cl);});        
            if (it != last_step_clusters.end()) {
                // This cluster was the same as it was last iteration (and it didn't advance) so just keep the same step size recommendation.
                cl.step_size = it->step_size;
                std::lock_guard guard(fine_min_final_time_lock);
                fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + it->step_size);
                printf("%i skipping time step selection because time step has already been picked for this cluster since last modification\n", cluster_i);
                return;
            }

            __itt_task_begin(domain, __itt_null, __itt_null, select_time_step_task);
            double step = _local_pde.selectTimeStep(cl);
            __itt_task_end(domain);
            {
                std::lock_guard guard(fine_min_final_time_lock);
                fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + step);
            }

            printf("Selected time: %f for cluster %i\n", step, cluster_i);

            for (Particle* p : cl.particles) {
                p->cluster_id = cluster_i;
            }
        });
    }
    task_group.wait();

    for (Particle* p : particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    // Splitting into two loops makes this point act like a sync point.
    //    Advantage - Only clusters that can advance are advanced.
    //    Disadvantage - We have to wait for all time step selections to be done.

    tbb::task_group step_task_group;

    // printf("Advance\n");
    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        step_task_group.run([&, cluster_i] {
            if (clusters[cluster_i].min_current_time > fine_min_final_time)
            {
                // printf("%i not advancing because another cluster too far behind\n", cluster_i);
                // This cluster can't advance because there is another cluster that is too far behind
                return;
            }

            if (clusters[cluster_i].is_static)
            {
                // printf("%i not advancing because it's static\n", cluster_i);
                // This cluster shouldn't be advanced as it's asleep or static
                return;
            }

            if (clusters[cluster_i].sleeping) {
                _local_pde.stepSleeping(clusters[cluster_i]);

                // printf("%i (sleeping) stepped forward %f\n", clusters[cluster_i].step_size, cluster_i);
            }
            else {
                // printf("Advancing from %f to %f (%f)\n", clusters[cluster_i].min_current_time, clusters[cluster_i].min_current_time + clusters[cluster_i].step_size, clusters[cluster_i].step_size);
                stepRecursive(clusters[cluster_i], true);
                for (Particle* p : clusters[cluster_i].particles) {
                    if (!p->is_static) {
                        // printf("Post step %i, last: %f, %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
                        break;
                    }
                }

                // printf("%i stepped forward %f\n", cluster_i, clusters[cluster_i].step_size);
            }
        });
    }

    step_task_group.wait();

    double min_current = 100000.0;
    for (Particle* p : particles) {
        if (!p->is_static) {
            assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
            min_current = std::min(min_current, p->current_state.getTime());
            // printf("Particle %i last: %f, current: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
        }
    }
    for (Particle* p : particles) {
        if (!p->is_static) {
           assert(p->last_state.getTime() <= min_current);    
        }
    }

    last_step_clusters.swap(clusters);

    __itt_task_end(domain);
}

#include "broad_phase_embree_cluster.h"

#include "../../model/particle.h"
#include "../../collision_detection/broad_phase_embree.h"
#include "../../collision_detection/separate_clusters.h"
#include "../../common/viewer.h"

using namespace Delta2;
using namespace strategy;

BroadPhaseEmbreeCluster::BroadPhaseEmbreeCluster(PDEStrategy& local_pde, common::Options& opt) :
                                                 BroadPhaseStrategy(local_pde, opt) {
    
}

void BroadPhaseEmbreeCluster::stepRecursive(Delta2::collision::Cluster& cluster, bool first_call) {
    printf("Step size: %f\n", cluster.step_size);
    double start_step = cluster.step_size;
    bool success = false;
    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         printf("stepRecursive  pre: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }
    success = _local_pde.step(cluster);
    // for (Particle* p : cluster.particles) {
    //     if (!p->is_static) {
    //         printf("stepRecursive post: %i, last time: %f, current time: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
    //         break;
    //     }
    // }
    if (!success) {
        printf("Fail\n");
        // Interpolate current_new back to 50% between last_start and current_start
        // Set future_new to current_start
        // Step forward from there to the current_start (recursive call)
        // Store current_new (to use in last_final)
        // Set timestep size to step_size_start / 2.0 and project new future state
        // Step forward (recursive call)
        // Step forward (recursive call)
        // Set last_final to the stored current_new and each of the particles last_time_step_size to min(last_time_step_size, step_size_start / 2.0)
        
        double new_min_time = 0.0;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                new_min_time = std::max(new_min_time, (p->last_state.getTime() + p->current_state.getTime()) / 2.0);
            }
        }

        cluster.min_current_time = new_min_time; // in the next stepRecursive call all the particles will get rolled back to this time
        printf("Rolling back to %f\n", cluster.min_current_time);

        double step_size;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                p->future_state = p->current_state;
                p->current_state = p->current_state.interpolate(p->last_state, new_min_time);
                // cluster.min_current_time = std::min(cluster.min_current_time, p->current_state.getTime());
                step_size = p->future_state.getTime() - new_min_time;
            }
        }

        cluster.step_size = step_size;
        stepRecursive(cluster, false);

        std::vector<State> last_state_final;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                last_state_final.push_back(p->current_state);
            }
        }

        cluster.step_size = start_step / 2.0;

        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                cluster.min_current_time = p->current_state.getTime();
                p->projectFutureState(cluster.step_size);
            }
        }

        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                printf("Step 2 %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
            }
        }
        stepRecursive(cluster, false);
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                cluster.min_current_time = p->current_state.getTime();
                printf("Step 3 %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
            }
        }
        stepRecursive(cluster, false);
        
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                printf("Post steps %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
            }
        }

        int i = 0;
        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                p->last_state = last_state_final[i++];
                p->last_time_step_size = std::min(p->last_time_step_size, start_step / 2.0);
            }
        }

        for (Particle* p : cluster.particles) {
            if (!p->is_static) {
                printf("End %i, last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
            }
        }
        
    }
    else {
        printf("Success\n");
    }

    cluster.step_size = start_step;
}

void BroadPhaseEmbreeCluster::step(model::ParticleHandler& particles) {
    for (Particle* p : particles)
    {
        p->last_time_step_size *= 2;
        p->last_time_step_size = std::min(p->last_time_step_size, (double) _opt.time_step_size);
        p->projectFutureState(p->last_time_step_size);
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    collision::BroadPhaseCollisions B;
    B = collision::broadPhaseEmbree(particles);

    std::vector<collision::Cluster> clusters = collision::separateCollisionClusters(B, particles);
    
    double min_final_time = std::numeric_limits<double>::infinity();

    for (collision::Cluster& cluster : clusters) {
        if (!cluster.is_static) {
            double cluster_min_time = cluster.min_current_time + cluster.step_size;
            min_final_time = std::min(min_final_time, cluster_min_time);
        }
    }

    // printf("Clusters: %i\n", clusters.size());

    double fine_min_final_time = std::numeric_limits<double>::infinity();

    // printf("Select time step sizes\n");
    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        collision::Cluster& cl = clusters[cluster_i];

        if (cl.min_current_time > min_final_time)
        {
            printf("%i skipping time step selection because min current time ahead of min final time\n", cluster_i);
            // This cluster can't advance because there is another cluster that is too far behind
            continue;
        }

        if (cl.is_static)
        {
            // printf("%i skipping time step selection because static cluster\n", cluster_i);
            // This cluster shouldn't be advanced as it's asleep or static
            continue;
        }

        if (cl.sleeping) {
            printf("%i skipping time step selection because sleeping cluster\n", cluster_i);
            fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + _opt.time_step_size);
            continue;
        }

        // TODO can this test ever pass?  Does the first check for min_current_time > min_final_time catch this?
        auto it = std::find_if(last_step_clusters.begin(), last_step_clusters.end(), [&cl](const collision::Cluster& last) {return !last.hasAdvanced(cl);});        
        if (it != last_step_clusters.end()) {
            // This cluster was the same as it was last iteration (and it didn't advance) so just keep the same step size recommendation.
            cl.step_size = it->step_size;
            fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + it->step_size);
            printf("%i skipping time step selection because time step has already been picked for this cluster since last modification\n", cluster_i);
            continue;
        }

        double step = _local_pde.selectTimeStep(cl);
        fine_min_final_time = std::min(fine_min_final_time, cl.min_current_time + step);

        printf("Selected time: %f for cluster %i\n", step, cluster_i);

        for (Particle* p : cl.particles) {
            p->cluster_id = cluster_i;
        }
    }

    for (Particle* p : particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    // Splitting into two loops makes this point act like a sync point.
    //    Advantage - Only clusters that can advance are advanced.
    //    Disadvantage - We have to wait for all time step selections to be done.

    // printf("Advance\n");
    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        if (clusters[cluster_i].min_current_time > fine_min_final_time)
        {
            // printf("%i not advancing because another cluster too far behind\n", cluster_i);
            // This cluster can't advance because there is another cluster that is too far behind
            continue;
        }

        if (clusters[cluster_i].is_static)
        {
            // printf("%i not advancing because it's static\n", cluster_i);
            // This cluster shouldn't be advanced as it's asleep or static
            continue;
        }

        if (clusters[cluster_i].sleeping) {
            _local_pde.stepSleeping(clusters[cluster_i]);

            // printf("%i (sleeping) stepped forward %f\n", clusters[cluster_i].step_size, cluster_i);
        }
        else {
            printf("Advancing from %f to %f (%f)\n", clusters[cluster_i].min_current_time, clusters[cluster_i].min_current_time + clusters[cluster_i].step_size, clusters[cluster_i].step_size);
            stepRecursive(clusters[cluster_i], true);
            for (Particle* p : clusters[cluster_i].particles) {
                if (!p->is_static) {
                    printf("Post step %i, last: %f, %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
                    break;
                }
            }

            // printf("%i stepped forward %f\n", cluster_i, clusters[cluster_i].step_size);
        }
    }

    double min_current = 100000.0;
    for (Particle* p : particles) {
        if (!p->is_static) {
            assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
            min_current = std::min(min_current, p->current_state.getTime());
            printf("Particle %i last: %f, current: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime());
        }
    }
    for (Particle* p : particles) {
        if (!p->is_static) {
           assert(p->last_state.getTime() <= min_current);    
        }
    }

    last_step_clusters.swap(clusters);
}

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

void BroadPhaseEmbreeCluster::step(model::ParticleHandler& particles) {
    for (Particle* p : particles)
    {
        p->last_time_step_size *= 2;
        p->last_time_step_size = std::min(p->last_time_step_size, (double) _opt.time_step_size);
        p->projectFutureState(p->last_time_step_size);
        assert(p->last_state.getTime() < p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
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

    printf("Clusters: %i\n", clusters.size());

    double fine_min_final_time = std::numeric_limits<double>::infinity();

    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        if (clusters[cluster_i].min_current_time > min_final_time)
        {
            // This cluster can't advance because there is another cluster that is too far behind
            continue;
        }

        if (clusters[cluster_i].is_static)
        {
            // This cluster shouldn't be advanced as it's asleep or static
            continue;
        }

        if (clusters[cluster_i].sleeping) {
            fine_min_final_time = std::min(fine_min_final_time, clusters[cluster_i].min_current_time + _opt.time_step_size);
        }

        double step = _local_pde.selectTimeStep(clusters[cluster_i]); 
        fine_min_final_time = std::min(fine_min_final_time, clusters[cluster_i].min_current_time + step);

        printf("Selected time: %f for cluster %i\n", step, cluster_i);
    }

    for (Particle* p : particles)
    {
        assert(p->last_state.getTime() < p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    // Splitting into two loops makes this point act like a sync point.
    //    Advantage - Only clusters that can advance are advanced.
    //    Disadvantage - We have to wait for all time step selections to be done.

    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++)
    {
        if (clusters[cluster_i].min_current_time > fine_min_final_time)
        {
            // This cluster can't advance because there is another cluster that is too far behind
            continue;
        }

        if (clusters[cluster_i].is_static)
        {
            // This cluster shouldn't be advanced as it's asleep or static
            continue;
        }

        if (clusters[cluster_i].sleeping) {
            _local_pde.stepSleeping(clusters[cluster_i]);

            printf("Sleeping cluster %i stepped forward %f\n", clusters[cluster_i].step_size, cluster_i);
        }
        else {
            _local_pde.step(clusters[cluster_i]);

            printf("Cluster %i stepped forward %f\n", clusters[cluster_i].step_size, cluster_i);
        }
    }
}

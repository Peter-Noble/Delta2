#pragma once

#include "broad_phase.h"
#include "../model/particle.h"
#include "../model/particle_handler.h"
#include "contact.h"
#include "contact_state.h"
#include "cluster.h"

namespace Delta2 {
    namespace collision {
        void separateContinuousCollisionClusters(BroadPhaseCollisions& broad_phase, std::vector<Delta2::Particle>& particles, double max_time, std::vector<std::vector<Particle*>>& cluster_particles_out, std::vector<BroadPhaseCollisions>& cluster_interactions_out, std::vector<bool>& sleeping);
        void separateCollisionClustersWithTimeStepSelection(BroadPhaseCollisions& broad_phase, std::vector<Delta2::Particle>& particles, double max_time, std::vector<std::vector<Particle*>>& cluster_particles_out, std::vector<BroadPhaseCollisions>& cluster_interactions_out, std::vector<double>& cluster_step_size_out, std::vector<bool>& sleeping);
        std::vector<Cluster> separateCollisionClusters(BroadPhaseCollisions& broad_phase, model::ParticleHandler& particles);
        void fineCollisionClustersWithTimeStepSelection(Cluster& cluster);
        void fineWitnessCollisionClustersWithTimeStepSelection(std::vector<Delta2::Particle>& particles, std::vector<std::vector<Delta2::Particle*>>& cluster_particles_out, std::vector<BroadPhaseCollisions>& cluster_interactions_out, std::vector<double>& cluster_step_size_out, std::vector<bool>& sleeping, std::vector<double>& min_current_time, collision::ContactStateCache& cache);
        std::vector<Cluster> separateClusterByTimestep(Cluster& cluster);
    }
}

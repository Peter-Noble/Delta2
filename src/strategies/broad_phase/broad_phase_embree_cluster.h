#pragma once

#include "broad_phase_strategy.h"

namespace Delta2 {
    namespace strategy {
        class BroadPhaseEmbreeCluster : public BroadPhaseStrategy {
        public:
            BroadPhaseEmbreeCluster(PDEStrategy& local_pde, common::Options& opt);
            void step(model::ParticleHandler& particles) override;
            void stepRecursive(Delta2::collision::Cluster& cluster, bool first_call, int depth=0);
        private:
            std::vector<collision::Cluster> last_step_clusters;
        };
    }
}

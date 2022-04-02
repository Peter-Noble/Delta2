#pragma once

#include "timestep_scheme.h"

namespace Delta2 {
    namespace timestepping {
        class ExplicitScheme : public TimestepScheme {
        public:
            ExplicitScheme(std::vector<Delta2::Particle>* ps, common::Options& o);
            // void step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& view_draws) override;
            void stepSetup() override;
            void broadPhase() override;
            void narrowPhase() override;
            void solveContacts() override;
            void integrate(uint32_t p_i, bool apply) override;
            void friction() override;
        protected:
            collision::BroadPhaseCollisions _broad_phase_collisions;
            std::vector<Eigen::Vector3d> _forces;
            std::vector<Eigen::Vector3d> _torques;
            std::vector<uint32_t> _counts;
        private:
            std::unordered_map<int, int> _global_to_local_ids;
        };
    }
}
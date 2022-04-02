#pragma once

#include "timestep_scheme.h"
#include "../collision_detection/contact_state.h"

namespace Delta2 {
    namespace timestepping {
        class ExplicitAdaptiveClustersWithWitnessesScheme : public TimestepScheme {
            public:
                ExplicitAdaptiveClustersWithWitnessesScheme(std::vector<Delta2::Particle>* ps, std::function<Eigen::Vector3d(const Particle&)> external_force, common::Options& o);
                void step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& view_draws) override;
                virtual ~ExplicitAdaptiveClustersWithWitnessesScheme() = default;
            private:
                double _last_time_step_size;
                std::vector<double> _sleep_candidates;
                collision::ContactStateCache _cache;
        };
    }
}

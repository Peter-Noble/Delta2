#pragma once

#include "timestep_scheme.h"

namespace Delta2 {
    namespace timestepping {
        class ImplicitScheme : public TimestepScheme {
            public:
                ImplicitScheme(std::vector<Delta2::Particle>* ps, std::function<Eigen::Vector3d(const Particle&)> external_force, common::Options& o);
                void step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& view_draws) override;
            private:

        };
    }
}

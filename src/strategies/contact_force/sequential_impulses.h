#pragma once

#include "contact_force_strategy.h"

namespace Delta2 {
    namespace strategy {
        class SequentialImpulses : public ContactForceStrategy {
        public:
            SequentialImpulses(FrictionStrategy& friction, common::Options& opt);
            bool solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) override;
        };
    }
}

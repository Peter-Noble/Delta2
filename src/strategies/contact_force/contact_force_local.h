#pragma once

#include "contact_force_strategy.h"

namespace Delta2 {
    namespace strategy {
        class ContactForceLocal : public ContactForceStrategy {
        public:
            ContactForceLocal(FrictionStrategy& friction, common::Options& opt);
            bool solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits, bool allow_fail=true) override;
        };
    }
}

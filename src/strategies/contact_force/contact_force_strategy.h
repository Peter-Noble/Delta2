#pragma once

#include "../friction/friction_strategy.h"
#include "../../collision_detection/cluster.h"
#include "../../collision_detection/contact.h"

namespace Delta2 {
    namespace strategy {
        class ContactForceStrategy {
        public:
            ContactForceStrategy(FrictionStrategy& friction, common::Options& opt);
            virtual bool solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) = 0;
            std::function<Eigen::Vector3d(const Particle &)> external_force;
        protected:
            FrictionStrategy& _friction;
            common::Options _opt;
        };
    }
}

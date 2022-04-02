#pragma once

#include "../../model/particle_handler.h"
#include "../../common/cli_options.h"
#include "../../collision_detection/cluster.h"

namespace Delta2 {
    namespace strategy {
        class PDEStrategy {
        public:
            PDEStrategy(common::Options& opt);
            // get the time step size, get contacts, solve contacts/friction, 
            virtual double selectTimeStep(collision::Cluster& cluster) = 0;
            virtual void step(collision::Cluster& cluster) = 0;
            void stepSleeping(collision::Cluster& cluster);
        protected:
            common::Options _opt;
        };
    }
}

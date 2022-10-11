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
            virtual bool step(collision::Cluster& cluster, bool allow_fail=true) = 0;
            void stepSleeping(collision::Cluster& cluster);
            virtual void printType() {
                printf("PDEStrategy");
            }
        protected:
            common::Options _opt;
        };
    }
}

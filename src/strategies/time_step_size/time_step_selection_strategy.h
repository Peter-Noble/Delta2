#pragma once

#include "../../collision_detection/cluster.h"
#include "../../common/cli_options.h"
#include "../../model/particle_handler.h"

namespace Delta2 {
    namespace strategy {
        class TimeStepSelectionStrategy {
        public:
            TimeStepSelectionStrategy(common::Options& opt);
            // Returns the desired time step size.  -1 if all static.  final_time is the expected time stamp after this step.
            virtual double selectTimeStep(collision::Cluster& cluster) = 0;
            virtual void init(model::ParticleHandler& ps);
        protected:
            common::Options _opt;
        };
    }
}

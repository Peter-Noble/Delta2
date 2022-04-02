#pragma once

#include "../PDE/PDE_strategy.h"
#include "../time_step_size/time_step_selection_strategy.h"
#include "../../common/cli_options.h"

namespace Delta2 {
    namespace strategy {
        class BroadPhaseStrategy {
        public:
            BroadPhaseStrategy(PDEStrategy& local_pde, common::Options& opt);
            virtual void step(model::ParticleHandler& particles) = 0;
        protected:
            PDEStrategy& _local_pde;
            common::Options _opt;
        };
    }
}

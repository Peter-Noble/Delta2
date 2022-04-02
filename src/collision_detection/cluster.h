#pragma once

#include "../model/particle_handler.h"
#include "broad_phase.h"


namespace Delta2 {
    namespace collision {
        struct Cluster {
            model::ParticleHandler particles;
            collision::BroadPhaseCollisions interations;
            double step_size;
            bool sleeping;
            double min_current_time;
            bool is_static;
        };
    }
}
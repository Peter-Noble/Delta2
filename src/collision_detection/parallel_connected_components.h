#pragma once

#include "colour_hits.h"
#include "broad_phase.h"
#include "../model/particle.h"
#include "../model/particle_handler.h"
#include "contact.h"
#include "contact_state.h"
#include "cluster.h"

#include "tbb/tbb.h"

namespace Delta2 {
    namespace collision {
        std::vector<uint32_t> simpleParallelConnectedComponents(collision::BroadPhaseCollisions& broad_phase, model::ParticleHandler& particles);
    }
}

#pragma once

#include <embree3/rtcore.h>
#include <stdio.h>
#include <limits>
#include <mutex>
#include <vector>

#include "../model/particle.h"
#include "broad_phase.h"
#include "../model/particle_handler.h"

#include "../embree_common/math/bbox.h"
#include "../embree_common/math/vec3fa.h"


namespace Delta2 {
    namespace collision {
        BroadPhaseCollisions broadPhaseEmbree(model::ParticleHandler& particles, double t);
        BroadPhaseCollisions broadPhaseEmbree(model::ParticleHandler& particles);
    }
}

#pragma once

#include "../common/primitive_geo.h"
#include "../common/basic_utils.h"
#include "../globals.h"
#include "../model/particle.h"
#include "../strategies/contact_force/contact_force_strategy.h"

namespace Delta2 {
    namespace scenarios {
        void make_scenario(std::vector<Delta2::Particle>& particles, Delta2::strategy::ContactForceStrategy& force_strategy);
    }
}

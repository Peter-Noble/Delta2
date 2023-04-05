#include "scenarios.h"

#include "../scenarios/zero_g_bundle.h"
#include "../scenarios/towers.h"
#include "../scenarios/start_intersecting.h"
#include "../scenarios/pairs.h"
#include "../scenarios/waterfall.h"
#include "../scenarios/fast_moving.h"
#include "../scenarios/sliding.h"
#include "../scenarios/tiny_tower.h"
#include "../scenarios/angle_sliding.h"
#include "../scenarios/big_pile.h"

void Delta2::scenarios::make_scenario(std::vector<Delta2::Particle>& particles, Delta2::strategy::ContactForceStrategy& force_strategy) {
    switch (globals::opt.scenario) {
    case 0:
        {
            scenarios::zero_g_bundle_orig(particles, force_strategy);
            break;
        }
    case 1:
        {
            scenarios::towers(particles, force_strategy);
            break;
        }
    case 2:
        {
            scenarios::start_intersecting(particles, force_strategy);
            break;
        }
    case 3:
        {
            scenarios::pairs(particles, force_strategy);
            break;
        }
    case 4:
        {
            scenarios::zero_g_bundle(particles, force_strategy);
            break;
        }
    case 5:
        {
            scenarios::waterfall(particles, force_strategy);
            break;
        }
    case 6:
        {
            scenarios::fast_moving(particles, force_strategy);
            break;
        }
    case 7:
        {
            scenarios::sliding(particles, force_strategy);
            break;
        }
    case 8:
        {
            scenarios::tiny_tower(particles, force_strategy);
            break;
        }
    case 9:
        {
            scenarios::angle_sliding(particles, force_strategy);
            break;
        }
    case 10:
        {
            scenarios::angle_sliding_first(particles, force_strategy);
            break;
        }
    case 11:
        {
            scenarios::pyramid(particles, force_strategy);
            break;
        }
    default:
        throw std::runtime_error("No scenario given");
    }
}

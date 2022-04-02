#include "PDE_strategy.h"
#include "../../model/integrate.h"

using namespace Delta2;
using namespace strategy;

PDEStrategy::PDEStrategy(common::Options& opt) :
                         _opt(opt) {
}

void PDEStrategy::stepSleeping(collision::Cluster& cluster) {
    for (Particle* p : cluster.particles) {
        model::staticApply(*p, cluster.step_size, true);
    }
}

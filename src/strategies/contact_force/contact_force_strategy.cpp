#include "contact_force_strategy.h"

using namespace Delta2;
using namespace strategy;

ContactForceStrategy::ContactForceStrategy(FrictionStrategy& friction,
                                           common::Options& opt) :
                                           _opt(opt),
                                           _friction(friction) {
    external_force = [](const Delta2::Particle& p) {
        // Eigen::Vector3d f = p.getMass() * Eigen::Vector3d({0, 0, -9.8});
        // return f; 
        return Eigen::Vector3d({0, 0, 0});
    };
}

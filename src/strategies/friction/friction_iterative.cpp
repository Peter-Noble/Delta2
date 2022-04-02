#include "friction_iterative.h"

#include "../../model/forces.h"

using namespace Delta2;
using namespace strategy;

FrictionIterative::FrictionIterative(common::Options& opt) :
                                     FrictionStrategy(opt) {

}

void FrictionIterative::solve(model::ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size) {
    Delta2::model::friction_solve(particles, hits, forces, torques, counts, external_force, step_size);
}
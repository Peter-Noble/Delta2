#include "friction_local.h"

using namespace Delta2;
using namespace strategy;

FrictionLocal::FrictionLocal(common::Options& opt) :
                             FrictionStrategy(opt) {

}

void FrictionLocal::solve(model::ParticleHandler& particles, std::vector<collision::Contact<double>>& hits, std::vector<Eigen::Vector3d>& forces, std::vector<Eigen::Vector3d>& torques, std::vector<int>& counts, std::function<Eigen::Vector3d(const Delta2::Particle&)> external_force, double step_size) {
    throw;
}

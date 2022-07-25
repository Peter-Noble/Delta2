#include "../collision_detection/contact.h"
#include "particle_handler.h"
#include <Eigen/Dense>

namespace Delta2 {
    namespace model {
        class FrictionSolver {
        public:
            FrictionSolver();
            FrictionSolver(ParticleHandler& particles,
                           std::vector<collision::Contact<double>>& input_hits);
            void solve(int add_iterations,
                       std::vector<Eigen::Vector3d>& forces,
                       std::vector<Eigen::Vector3d>& torques,
                       std::vector<State>& future_states,
                       double step_size);
            std::vector<Eigen::Vector3d> friction_forces;
            std::vector<Eigen::Vector3d> friction_torques;
        private:
            ParticleHandler* ph;
            std::vector<collision::Contact<double>>* hits;
            std::vector<Eigen::Vector3d> start_rel_vels;
            std::vector<double> max_force;
            std::vector<double> friction_factor;
            std::vector<Eigen::Vector3d> hit_friction;

            double adjust_by = 0.1;
            double rate = 0.95;
            int iterations = 0;
        };
    }
}

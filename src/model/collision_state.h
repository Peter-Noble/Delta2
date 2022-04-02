#pragma once

#include <vector>
#include "particle.h"

namespace Delta2 {
    namespace collision {
		struct NodeState {
			int node_id;
			bool expand;
			bool reduce;
			bool newly_expanded;
		};
		struct ActiveSet {
			std::vector<NodeState> states;
		};

		class PairState {
		public:
			PairState(Particle* A, Particle* B);
			void updateActiveSets();
		private:
			Particle *particle_A, *particle_B;
			ActiveSet A, B;

			Eigen::Vector3d a_last_force;
			Eigen::Vector3d a_avg_force;
			Eigen::Vector3d a_last_update_force;
			Eigen::Vector3d a_last_torque;
			Eigen::Vector3d a_avg_torque;
			Eigen::Vector3d a_last_update_torque;
			Eigen::Vector3d a_force;
			Eigen::Vector3d a_torque;

			Eigen::Vector3d b_last_force;
			Eigen::Vector3d b_avg_force;
			Eigen::Vector3d b_last_update_force;
			Eigen::Vector3d b_last_torque;
			Eigen::Vector3d b_avg_torque;
			Eigen::Vector3d b_last_update_torque;
			Eigen::Vector3d b_force;
			Eigen::Vector3d b_torque;
		};

		struct DeferredCompare {
			DeferredCompare(const model::Node&, const model::Node&);
			const model::Node &a;
			const model::Node &b;
		};

		struct Island {
			std::vector<PairState> pairs;
		};
	}
}

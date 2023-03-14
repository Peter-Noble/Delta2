#pragma once

#include <Eigen/Dense>

#include "../common/triangle.h"
#include "../common/aligned.h"
#include "contact.h"
#include "../model/surrogate_components.h"
#include "../model/collision_state.h"
#include "../model/surrogate_tree.h"
#include "../model/particle.h"

// TODO this is currently just a renamed copy of compare.  Use the sympy_bucket_deferred from hybrid_dist.cpp

namespace Delta2 {
	namespace collision {
        std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> _hybridDists(
			std::vector<Delta2::common::Triangle<float>>& As,
			std::vector<Delta2::common::Triangle<float>>& Bs,
			std::vector<float>& max_error);

		void findContactsBucketHybridDeferred(
            const std::vector<DeferredCompare>& bucket_pairs,
            Particle& particle_a,
            Particle& particle_b,
            Eigen::Matrix4d a_trans,
            Eigen::Matrix4d b_trans,
            std::vector<Contact<double>>& hits,
			std::vector<int>& pair_used_out);
	}
}
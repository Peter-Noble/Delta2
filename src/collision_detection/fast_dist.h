#pragma once

#include <Eigen/Dense>

#include "../common/triangle.h"
#include "contact.h"
#include "../model/surrogate_components.h"
#include "../model/collision_state.h"
#include "../model/surrogate_tree.h"
#include "../model/particle.h"

namespace Delta2 {
	namespace collision {
		std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> _fastInnerDists(
			Delta2::common::Triangle<float>& a_tri,
			std::vector<Delta2::common::Triangle<float>>& b_tris);

		template<int branching, typename real>
		void findContactsBucketFast(std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits) {
			float search_dist = a.geo_eps + b.geo_eps;

			for (DeferredCompare& dc : bucket_pairs) {
				std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
				std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

				std::vector<real> a_epss = a_bucket->getEps();
				std::vector<real> b_epss = b_bucket->getEps();
				std::vector<real> a_inner_epss = a_bucket->getInnerEps();
				std::vector<real> b_inner_epss = b_bucket->getInnerEps();

				int a_i = 0;
				for (common::Triangle<real>& a_tri : a_bucket->getTriangles()) {
					common::Triangle<float> A = a_tri.template cast<float>().transformed(a.future_state.getTransformation().cast<float>());

					std::vector<Delta2::common::Triangle<float>> Bs;
					for(common::Triangle<real> b_tri : b_bucket->getTriangles()) {
						Bs.push_back(b_tri.template cast<float>().transformed(b.future_state.getTransformation().cast<float>()));
					}

					int b_i = 0;
					std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = _fastInnerDists(A, Bs);
					for (auto& [dist, P, Q] : result) {
						if (dist < search_dist + a_epss[a_i] + b_epss[b_i]) {
							Contact<real> ci(P.cast<real>(), Q.cast<real>(), a.geo_eps, b.geo_eps, 0.0, 0.0, a, b);
							hits.push_back(ci);
						}
						b_i++;
					}
					a_i++;
				}
			}
		};
	}
}

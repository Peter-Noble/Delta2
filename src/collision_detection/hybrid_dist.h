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
        template<int branching, class real>
		void findContactsBucketHybrid(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
			float search_dist = a.geo_eps + b.geo_eps;
            pair_used_out.reserve(bucket_pairs.size());
            pair_used_out.clear();

			for (const DeferredCompare& dc : bucket_pairs) {
				std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
				std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

				std::vector<real> a_epss = a_bucket->getEps();
				std::vector<real> b_epss = b_bucket->getEps();
				std::vector<real> a_inner_epss = a_bucket->getInnerEps();
				std::vector<real> b_inner_epss = b_bucket->getInnerEps();

                // TODO loop over vertices and edges rather than triangles to reduce duplicated checks
                std::vector<common::Triangle<real>> a_tris = a_bucket->getTriangles();
                std::vector<common::Triangle<real>> b_tris = b_bucket->getTriangles();
                for (int a_i = 0; a_i < a_tris.size(); a_i++) {
                    common::Triangle<real>& a_tri = a_tris[a_i];
					common::Triangle<real> A = a_tri.transformed(a.future_state.getTransformation().cast<real>());

                    for (int b_i = 0; b_i < b_tris.size(); b_i++) {
                        common::Triangle<real>& b_tri = b_tris[b_i];
                        common::Triangle<real> B = b_tri.transformed(b.future_state.getTransformation().cast<real>());

                        auto [dist, P, Q] = _hybridDist(A, B);
						if (dist < search_dist + a_epss[a_i] + b_epss[b_i]) {
							Contact<real> ci(P, Q, a.geo_eps + a_epss[a_i], b.geo_eps + b_epss[b_i], a_inner_epss[a_i], b_inner_epss[b_i], a, b);
							hits.push_back(ci);
						}
                    }
				}
                pair_used_out.push_back(hits.size());
			}
		};

        template<int branching, class real>
		void findContactsBucketHybridCurrent(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
			float search_dist = a.geo_eps + b.geo_eps;
            pair_used_out.reserve(bucket_pairs.size());
            pair_used_out.clear();

			for (const DeferredCompare& dc : bucket_pairs) {
				std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
				std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

				std::vector<real> a_epss = a_bucket->getEps();
				std::vector<real> b_epss = b_bucket->getEps();
				std::vector<real> a_inner_epss = a_bucket->getInnerEps();
				std::vector<real> b_inner_epss = b_bucket->getInnerEps();

                // TODO loop over vertices and edges rather than triangles to reduce duplicated checks
                std::vector<common::Triangle<real>> a_tris = a_bucket->getTriangles();
                std::vector<common::Triangle<real>> b_tris = b_bucket->getTriangles();
                for (int a_i = 0; a_i < a_tris.size(); a_i++) {
                    common::Triangle<real>& a_tri = a_tris[a_i];
					common::Triangle<real> A = a_tri.transformed(a.current_state.getTransformation().cast<real>());

                    for (int b_i = 0; b_i < b_tris.size(); b_i++) {
                        common::Triangle<real>& b_tri = b_tris[b_i];
                        common::Triangle<real> B = b_tri.transformed(b.current_state.getTransformation().cast<real>());

                        auto [dist, P, Q] = _hybridDist(A, B);
						if (dist < search_dist + a_epss[a_i] + b_epss[b_i]) {
							Contact<real> ci(P, Q, a.geo_eps + a_epss[a_i], b.geo_eps + b_epss[b_i], a_inner_epss[a_i], b_inner_epss[b_i], a, b);
							hits.push_back(ci);
						}
                    }
				}
                pair_used_out.push_back(hits.size());
			}
		};

        template<int branching, class real>
		void findContactsBucketHybridLast(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
			float search_dist = a.geo_eps + b.geo_eps;
            pair_used_out.reserve(bucket_pairs.size());
            pair_used_out.clear();

			for (const DeferredCompare& dc : bucket_pairs) {
				std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
				std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

				std::vector<real> a_epss = a_bucket->getEps();
				std::vector<real> b_epss = b_bucket->getEps();
				std::vector<real> a_inner_epss = a_bucket->getInnerEps();
				std::vector<real> b_inner_epss = b_bucket->getInnerEps();

                // TODO loop over vertices and edges rather than triangles to reduce duplicated checks
                std::vector<common::Triangle<real>> a_tris = a_bucket->getTriangles();
                std::vector<common::Triangle<real>> b_tris = b_bucket->getTriangles();
                for (int a_i = 0; a_i < a_tris.size(); a_i++) {
                    common::Triangle<real>& a_tri = a_tris[a_i];
					common::Triangle<real> A = a_tri.transformed(a.last_state.getTransformation().cast<real>());

                    for (int b_i = 0; b_i < b_tris.size(); b_i++) {
                        common::Triangle<real>& b_tri = b_tris[b_i];
                        common::Triangle<real> B = b_tri.transformed(b.last_state.getTransformation().cast<real>());

                        auto [dist, P, Q] = _hybridDist(A, B);
						if (dist < search_dist + a_epss[a_i] + b_epss[b_i]) {
							Contact<real> ci(P, Q, a.geo_eps + a_epss[a_i], b.geo_eps + b_epss[b_i], a_inner_epss[a_i], b_inner_epss[b_i], a, b);
							hits.push_back(ci);
						}
                    }
				}
                pair_used_out.push_back(hits.size());
			}
		};
	}
}
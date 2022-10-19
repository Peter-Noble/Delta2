#pragma once

#include <Eigen/Dense>

#include "../common/triangle.h"
#include "contact.h"
#include "../model/collision_state.h"
#include "../common/viewer.h"

namespace Delta2 {
	namespace collision {
        template<typename real>
		std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>> _comparisonDist(Delta2::common::Triangle<real>& a_tri, Delta2::common::Triangle<real>& b_tri) {
            Eigen::Vector<real, 3> P_closest, Q_closest, P_tmp, Q_tmp;
            real dist_closest = std::numeric_limits<real>::infinity();
            real dist;
            bool isInTriangle;
            
            P_tmp = a_tri.projectToPlane(b_tri.A, isInTriangle);
            if (isInTriangle) {
                dist = (P_tmp - b_tri.A).norm();
                if (dist < dist_closest) {
                    P_closest = P_tmp;
                    Q_closest = b_tri.A;
                    dist_closest = dist;
                }
            }

            P_tmp = a_tri.projectToPlane(b_tri.B, isInTriangle);
            if (isInTriangle) {
                dist = (P_tmp - b_tri.B).norm();
                if (dist < dist_closest) {
                    P_closest = P_tmp;
                    Q_closest = b_tri.B;
                    dist_closest = dist;
                }
            }

            P_tmp = a_tri.projectToPlane(b_tri.C, isInTriangle);
            if (isInTriangle) {
                dist = (P_tmp - b_tri.C).norm();
                if (dist < dist_closest) {
                    P_closest = P_tmp;
                    Q_closest = b_tri.C;
                    dist_closest = dist;
                }
            }


            Q_tmp = b_tri.projectToPlane(a_tri.A, isInTriangle);
            if (isInTriangle) {
                dist = (a_tri.A - Q_tmp).norm();
                if (dist < dist_closest) {
                    P_closest = a_tri.A;
                    Q_closest = Q_tmp;
                    dist_closest = dist;
                }
            }

            Q_tmp = b_tri.projectToPlane(a_tri.B, isInTriangle);
            if (isInTriangle) {
                dist = (a_tri.B - Q_tmp).norm();
                if (dist < dist_closest) {
                    P_closest = a_tri.B;
                    Q_closest = Q_tmp;
                    dist_closest = dist;
                }
            }

            Q_tmp = b_tri.projectToPlane(a_tri.C, isInTriangle);
            if (isInTriangle) {
                dist = (a_tri.C - Q_tmp).norm();
                if (dist < dist_closest) {
                    P_closest = a_tri.C;
                    Q_closest = Q_tmp;
                    dist_closest = dist;
                }
            }

            dist = a_tri.AB().closest(b_tri.AB(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.AB().closest(b_tri.BC(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.AB().closest(b_tri.CA(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.BC().closest(b_tri.AB(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.BC().closest(b_tri.BC(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.BC().closest(b_tri.CA(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }
            dist = a_tri.CA().closest(b_tri.AB(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.CA().closest(b_tri.BC(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            dist = a_tri.CA().closest(b_tri.CA(), P_tmp, Q_tmp);
            if (dist < dist_closest) {
                P_closest = P_tmp;
                Q_closest = Q_tmp;
                dist_closest = dist;
            }

            real t; 
            bool isIntersecting;
            bool t_inf;
            a_tri.intersectSegment(b_tri.AB(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = b_tri.AB().lerp(t);
                Q_closest = b_tri.AB().lerp(t);
                dist_closest = 0.0;
            }

            a_tri.intersectSegment(b_tri.BC(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = b_tri.BC().lerp(t);
                Q_closest = b_tri.BC().lerp(t);
                dist_closest = 0.0;
            }

            a_tri.intersectSegment(b_tri.CA(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = b_tri.CA().lerp(t);
                Q_closest = b_tri.CA().lerp(t);
                dist_closest = 0.0;
            }


            b_tri.intersectSegment(a_tri.AB(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = a_tri.AB().lerp(t);
                Q_closest = a_tri.AB().lerp(t);
                dist_closest = 0.0;
            }

            b_tri.intersectSegment(a_tri.BC(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = a_tri.BC().lerp(t);
                Q_closest = a_tri.BC().lerp(t);
                dist_closest = 0.0;
            }

            b_tri.intersectSegment(a_tri.CA(), t, t_inf, isIntersecting);
            if (isIntersecting) {
                P_closest = a_tri.CA().lerp(t);
                Q_closest = a_tri.CA().lerp(t);
                dist_closest = 0.0;
            }

            return std::make_tuple((P_closest - Q_closest).norm(), P_closest, Q_closest);
        }
        template std::tuple<float, Eigen::Vector<float, 3>, Eigen::Vector<float, 3>> _comparisonDist(Delta2::common::Triangle<float>&, Delta2::common::Triangle<float>&);
        template std::tuple<double, Eigen::Vector<double, 3>, Eigen::Vector<double, 3>> _comparisonDist(Delta2::common::Triangle<double>&, Delta2::common::Triangle<double>&);

        template<typename real>
        std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> _comparisonInnerDists(Delta2::common::Triangle<real>& a_tri, std::vector<Delta2::common::Triangle<real>>& b_tris) {
            std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> result;
            result.reserve(b_tris.size());
            for (auto& b_tri : b_tris) {
                result.push_back(_comparisonDist(a_tri, b_tri));
            }
            return result;
        }
        template std::vector<std::tuple<float, Eigen::Vector<float, 3>, Eigen::Vector<float, 3>>> _comparisonInnerDists(Delta2::common::Triangle<float>&, std::vector<Delta2::common::Triangle<float>>&);
        template std::vector<std::tuple<double, Eigen::Vector<double, 3>, Eigen::Vector<double, 3>>> _comparisonInnerDists(Delta2::common::Triangle<double>&, std::vector<Delta2::common::Triangle<double>>&);

        template<int branching, class real>
		void findContactsBucketComparison(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
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

                        auto [dist, P, Q] = _comparisonDist(A, B);
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
		void findContactsBucketComparisonCurrent(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
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

                        auto [dist, P, Q] = _comparisonDist(A, B);
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
		void findContactsBucketComparisonLast(const std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<Contact<real>>& hits, std::vector<int>& pair_used_out) {
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

                        auto [dist, P, Q] = _comparisonDist(A, B);
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

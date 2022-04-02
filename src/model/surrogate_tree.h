#pragma once

#include <vector>
#include <memory>
#include <queue>
//#include "../common/common.h"
#include "../common/utils.h"
#include "surrogate_components.h"
#include "construct_surrogate_tree.h"

namespace Delta2 {
	namespace model {
		enum class SurrogateTreeConstructionMethod { KMEAN, QSLIM };

		class SurrogateTree {
		public:
			SurrogateTree();
			SurrogateTree(const Eigen::Matrix<double, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt, SurrogateTreeConstructionMethod cm=SurrogateTreeConstructionMethod::KMEAN);
			const Node& getNode(int id);
			std::shared_ptr<Bucket> getBucket(int id);
			int getNumNodes();
		private:
			Node& getNodeMutable(int id);
			void buildSurrogateTreeKMean(const Eigen::Matrix<double, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt);
			int addNode(Node N);
			int addBucket(std::shared_ptr<Bucket> B);
			std::vector<Node> _nodes;
			std::vector<std::shared_ptr<Bucket>> _buckets;
		};

		/*template<typename real, int branching, int bucket_size>
		class SurrogateTreeTemplate {
			public:
				SurrogateTreeTemplate() {

				};

				SurrogateTreeTemplate(const Eigen::Matrix<real, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt, SurrogateTreeConstructionMethod cm=SurrogateTreeConstructionMethod::KMEAN) {
					switch (cm) {
						case (SurrogateTreeConstructionMethod::KMEAN):
							buildSurrogateTreeKMean(V, F, opt);
							break;
						case (SurrogateTreeConstructionMethod::QSLIM):
							// TODO
							break;
					}
				};
				
				const Node& getNode(int id) {
					// Unlike the buckets, nodes are small so all nodes should be in memory for quick access
					return getNodeMutable(id);
				};
				
				std::shared_ptr<Bucket> getBucket(int id) {
					// For super large meshes this could involved loading this bucket off disk
					return _buckets[id];
				};

				int getNumNodes() {
					return _nodes.size();
				}

				// bool compareNode(const Node& a, const SurrogateTreeTemplate<real, branching, bucket_size>& other, const Node& b) const {
				// 	std::shared_ptr<model::Bucket> a_b = getBucket(a.bucket_id);
				// 	std::shared_ptr<model::Bucket> b_b = other.getBucket(b.bucket_id);

				// 	collision::findContactsBucketComparison()
				// }
			private:
				Node& getNodeMutable(int id) {
					// Unlike the buckets, nodes are small so all nodes should be in memory for quick access
					return _nodes[id];
				};

				void buildSurrogateTreeKMean(const Eigen::Matrix<real, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt) {
					std::vector<common::Triangle<real>> tris = common::toTriangles<real>(V, F, F.rows());
					std::vector<int> tri_ids;
					for (int i = 0; i < tris.size(); i++) {
						tri_ids.push_back(i);
					}
					std::queue<std::tuple<std::vector<int>, int, int>> in_progress; // grp, depth, parent
					in_progress.push(std::make_tuple(tri_ids, 0, -1));

					while (!in_progress.empty()) {
						auto [ids, depth, parent] = in_progress.front();
						in_progress.pop();
						int node_id;
						if (ids.size() < bucket_size) {
							std::vector<common::Triangle<real>> bucket_tris;
							for (int id : ids) {
								bucket_tris.push_back(tris[id]);
							}
							std::shared_ptr<BucketSoup> sp = std::make_shared<BucketSoup>(bucket_tris);
							int bucket_id = addBucket(std::static_pointer_cast<Bucket>(sp));

							model::Node n;
							n.bucket_id = bucket_id;
							n.is_surface = false;
							n.is_inner = false;
							n.depth = depth;
							n.parent_id = parent;
							n.child_id_first = -1;
							n.num_children = 0;

							node_id = addNode(n);
						} else {
							std::vector<common::Triangle<real>> surrogate_tris;
							std::vector<std::vector<int>> surrogate_grp_ids;
							divideAndFit(tris, ids, branching, true, opt, surrogate_tris, surrogate_grp_ids);

							std::vector<real> inner_epss;
							std::vector<real> epss;

							for (int i = 0; i < surrogate_grp_ids.size(); i++) {
								real inner_eps;
								real eps;
								Delta2::model::calcEpss(surrogate_tris[i], tris, surrogate_grp_ids[i], inner_eps, eps);
								inner_epss.push_back(inner_eps);
								epss.push_back(eps);
							}

							std::shared_ptr<BucketSoup> sp = std::make_shared<BucketSoup>(surrogate_tris);
							sp->setEps(epss);
							sp->setInnerEps(inner_epss);
							int bucket_id = addBucket(std::static_pointer_cast<Bucket>(sp));

							model::Node n;
							n.bucket_id = bucket_id;
							n.is_surface = false;
							n.is_inner = true;
							n.depth = depth;
							n.parent_id = parent;
							n.child_id_first = -1;  // We'll set this when processing the first node that's a child of this one
							n.num_children = surrogate_tris.size();

							node_id = addNode(n);

							for (int i = 0; i < surrogate_tris.size(); i++) {
								Delta2::common::Triangle<real> s_tri = surrogate_tris[i];
								std::vector<int> grp_ids = surrogate_grp_ids[i];
								in_progress.push(std::make_tuple(grp_ids, depth + 1, node_id));
							}
						}
						if (parent != -1 && getNode(parent).child_id_first == -1) {
							getNodeMutable(parent).child_id_first = node_id;
						}
					}
				};
				int addNode(Node N) {
					_nodes.push_back(N);
					_nodes.back().node_id = _nodes.size() - 1;
					return _nodes.back().node_id;
				}
				int addBucket(std::shared_ptr<Bucket> B) {
					_buckets.push_back(B);
					return _buckets.size() - 1;
				};
				std::vector<Node> _nodes;
				std::vector<std::shared_ptr<Bucket>> _buckets;
		};

		using SurrogateTree = SurrogateTreeTemplate<double, 8, 8>; */
	}
}


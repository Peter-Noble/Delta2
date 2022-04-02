#include "surrogate_tree.h"

namespace Delta2 {
    namespace model {
        SurrogateTree::SurrogateTree() {

        };

        SurrogateTree::SurrogateTree(const Eigen::Matrix<double, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt, SurrogateTreeConstructionMethod cm) {
            switch (cm) {
                case (SurrogateTreeConstructionMethod::KMEAN):
                    buildSurrogateTreeKMean(V, F, opt);
                    break;
                case (SurrogateTreeConstructionMethod::QSLIM):
                    // TODO
                    break;
            }
        };
        
        const Node& SurrogateTree::getNode(int id) {
            // Unlike the buckets, nodes are small so all nodes should be in memory for quick access
            return getNodeMutable(id);
        };
        
        std::shared_ptr<Bucket> SurrogateTree::getBucket(int id) {
            // For super large meshes this could involved loading this bucket off disk
            return _buckets[id];
        };

        int SurrogateTree::getNumNodes() {
            return _nodes.size();
        }

        Node& SurrogateTree::getNodeMutable(int id) {
            // Unlike the buckets, nodes are small so all nodes should be in memory for quick access
            return _nodes[id];
        };

        void SurrogateTree::buildSurrogateTreeKMean(const Eigen::Matrix<double, -1, -1>& V, const Eigen::MatrixXi& F, common::Options& opt) {
            std::vector<common::Triangle<double>> tris = common::toTriangles<double>(V, F, F.rows());
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
                if (ids.size() < 8) {
                    std::vector<common::Triangle<double>> bucket_tris;
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
                    std::vector<common::Triangle<double>> surrogate_tris;
                    std::vector<std::vector<int>> surrogate_grp_ids;
                    divideAndFit(tris, ids, 8, true, opt, surrogate_tris, surrogate_grp_ids);

                    std::vector<double> inner_epss;
                    std::vector<double> epss;

                    for (int i = 0; i < surrogate_grp_ids.size(); i++) {
                        double inner_eps;
                        double eps;
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
                        Delta2::common::Triangle<double> s_tri = surrogate_tris[i];
                        std::vector<int> grp_ids = surrogate_grp_ids[i];
                        in_progress.push(std::make_tuple(grp_ids, depth + 1, node_id));
                    }
                }
                if (parent != -1 && getNode(parent).child_id_first == -1) {
                    getNodeMutable(parent).child_id_first = node_id;
                }
            }
        };
        int SurrogateTree::addNode(Node N) {
            _nodes.push_back(N);
            _nodes.back().node_id = _nodes.size() - 1;
            return _nodes.back().node_id;
        }
        int SurrogateTree::addBucket(std::shared_ptr<Bucket> B) {
            _buckets.push_back(B);
            return _buckets.size() - 1;
        };
    }
}

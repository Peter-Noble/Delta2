#include "surrogate_components.h"
#include "../common/viewer.h"
#include "../common/convex_hull.h"

namespace Delta2 {
    namespace model {
        // BucketConnected
        BucketConnected::BucketConnected(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<int, -1, -1>& F) {
            if (V.rows() > 3*8) {
                throw std::invalid_argument("Bucket can't have this many vertices");
            }
            if (F.rows() > 8) {
                throw std::invalid_argument("Bucket can't have this many faces");
            }

            _vertices = V;
            _faces = F;
            _num_tris = F.rows();

            for (int i = 0; i < _num_tris; i++) {
                _eps[i] = 0.0;
                _inner_eps[i] = 0.0;
            }

            _hull = common::convex_hull(getVertices());
        };
        std::vector<common::Triangle<double>> BucketConnected::getTriangles() {
            return common::toTriangles(_vertices, _faces, _num_tris);
        };
        std::vector<common::Edge<double>> BucketConnected::getEdges() {
            Eigen::MatrixXi edges;
            igl::edges(_faces, edges);

            std::vector<common::Edge<double>> result;
            for (int e = 0; e < edges.rows(); e++) {
                Eigen::Vector<double, 3> A(_vertices(edges(e, 0), 0), _vertices(edges(e, 0), 1), _vertices(edges(e, 0), 2));
                Eigen::Vector<double, 3> B(_vertices(edges(e, 1), 0), _vertices(edges(e, 1), 1), _vertices(edges(e, 1), 2));
                result.push_back(common::Edge<double>(A, B));
            }
            return result;
        };
        std::vector<Eigen::Vector<double, 3>> BucketConnected::getVertices() {
            std::vector<Eigen::Vector<double, 3>> result;
            for (int v = 0; v < _vertices.rows(); v++) {
                Eigen::Vector<double, 3> V(_vertices(v, 0), _vertices(v, 1), _vertices(v, 2));
                result.push_back(V);
            }
            return result;
        }
        common::Triangle<double> BucketConnected::getTriangle(uint32_t n) {
            if (n < _num_tris) {
                Eigen::Vector<double, 3> A(_vertices(_faces(n,0),0), _vertices(_faces(n,0),1), _vertices(_faces(n,0),2));
                Eigen::Vector<double, 3> B(_vertices(_faces(n,1),0), _vertices(_faces(n,1),1), _vertices(_faces(n,1),2));
                Eigen::Vector<double, 3> C(_vertices(_faces(n,2),0), _vertices(_faces(n,2),1), _vertices(_faces(n,2),2));
                return common::Triangle(A, B, C);
            }
            throw std::invalid_argument("Not that many triangles in this bucket");
        };
        void BucketConnected::setEps(std::vector<double>& eps_in) {
            std::copy(eps_in.begin(), eps_in.end(), _eps.begin());
        }
        void BucketConnected::setInnerEps(std::vector<double>& inner_eps_in) {
            std::copy(inner_eps_in.begin(), inner_eps_in.end(), _inner_eps.begin());
        }
        std::vector<double> BucketConnected::getEps() {
            std::vector<double> result;
            for (int i = 0; i < _num_tris; i++) {
                result.push_back(_eps[i]);
            }
            return result;
        }
        std::vector<double> BucketConnected::getInnerEps() {
            std::vector<double> result;
            for (int i = 0; i < _num_tris; i++) {
                result.push_back(_inner_eps[i]);
            }
            return result;
        }
        void BucketConnected::printBucketConfig() {
            printf("Connected bucket\n");
        }
        std::vector<common::Triangle<double>> BucketConnected::getConvexHull() {
            return _hull;
        }
        BucketConfig BucketConnected::getConfig() {
            return BucketConfig::CONNECTED;
        };
        double BucketConnected::getMaxEps() {
            return *std::max_element(_eps.begin(), _eps.end());
        }
    
    
    
        // Bucket soup
        BucketSoup::BucketSoup(std::vector<common::Triangle<double>>& tris) {
            if (tris.size() > 8) {
                throw std::invalid_argument("A bucket can only hold so many triangles.");
            }
            for (int i = 0; i < std::min(8, (int)tris.size()); i++) {
                _triangles[i] = tris[i];
            }
            _num_tris = tris.size();

            for (int i = 0; i < _num_tris; i++) {
                _eps[i] = 0.0;
                _inner_eps[i] = 0.0;
            }

            _hull = common::convex_hull(getVertices());
        };
        BucketSoup::BucketSoup(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<int, -1, -1>& F) {
            std::vector<common::Triangle<double>> tris = common::toTriangles(V, F, F.rows());
            std::copy(tris.begin(), tris.end(), _triangles.begin());
            _num_tris = F.rows();

            for (int i = 0; i < _num_tris; i++) {
                _eps[i] = 0.0;
                _inner_eps[i] = 0.0;
            }
            
            _hull = common::convex_hull(getVertices());
        };
        std::vector<common::Triangle<double>> BucketSoup::getTriangles() {
            std::vector<common::Triangle<double>> r(_triangles.begin(), _triangles.begin() + _num_tris);
            return r;
        };
        std::vector<common::Edge<double>> BucketSoup::getEdges() {
            std::vector<common::Edge<double>> result;
            for (int i = 0; i < _num_tris; i++) {
                const common::Triangle<double>& T = _triangles[i];
                result.push_back(T.AB());
                result.push_back(T.BC());
                result.push_back(T.CA());
            }
            return result;
        };
        std::vector<Eigen::Vector<double, 3>> BucketSoup::getVertices() {
            std::vector<Eigen::Vector<double, 3>> result;
            for (int i = 0; i < _num_tris; i++) {
                for (int v = 0; v < 3; v++) {
                    result.push_back(_triangles[i][v]);
                }
            }
            return result;
        }
        common::Triangle<double> BucketSoup::getTriangle(uint32_t n) {
            if (n < _num_tris) {
                return _triangles[n];
            }
            throw std::invalid_argument("Not that many triangles in this bucket");
        };
        void BucketSoup::setEps(std::vector<double>& eps_in) {
            for (int i = 0; i < std::min(_num_tris, (int)eps_in.size()); i++) {
                _eps[i] = eps_in[i];
            }
        }
        void BucketSoup::setInnerEps(std::vector<double>& inner_eps_in) {
            for (int i = 0; i < std::min(_num_tris, (int)inner_eps_in.size()); i++) {
                _inner_eps[i] = inner_eps_in[i];
            }
        }
        std::vector<double> BucketSoup::getEps() {
            std::vector<double> result;
            for (int i = 0; i < _num_tris; i++) {
                result.push_back(_eps[i]);
            }
            return result;
        }
        std::vector<double> BucketSoup::getInnerEps() {
            std::vector<double> result;
            for (int i = 0; i < _num_tris; i++) {
                result.push_back(_inner_eps[i]);
            }
            return result;
        }
        void BucketSoup::printBucketConfig() {
            printf("Soup bucket\n");
        }
        std::vector<common::Triangle<double>> BucketSoup::getConvexHull() {
            throw; // TODO
        }
        BucketConfig BucketSoup::getConfig() {
            return BucketConfig::SOUP;
        }
        double BucketSoup::getMaxEps() {
            return *std::max_element(_eps.begin(), _eps.begin()+_num_tris);
        }
    }
}

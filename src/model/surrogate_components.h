#pragma once

#include <vector>
#include "../common/triangle.h"
#include "igl/edges.h"
#include "igl/copyleft/cgal/convex_hull.h"
#include "../common/utils.h"

namespace Delta2 {
    namespace model {
		struct Node {
			uint32_t node_id;
			uint32_t bucket_id;
			bool is_surface, is_inner;
			uint32_t depth;
			uint32_t parent_id;
			uint32_t child_id_first;
			uint32_t num_children;
			// TODO add min/max eps here so we can do something like: get active set for eps=0.01 (eg for a multigrid fluid sim)
		};

		enum BucketConfig { SOUP, CONNECTED };

		class Bucket {
			public:
				virtual ~Bucket() {};
				virtual std::vector<common::Triangle<double>> getTriangles() = 0;
				virtual std::vector<common::Edge<double>> getEdges() = 0;
				virtual std::vector<Eigen::Vector<double, 3>> getVertices() = 0;
				virtual common::Triangle<double> getTriangle(uint32_t n) = 0;
				virtual void setEps(std::vector<double>&) = 0;
				virtual void setInnerEps(std::vector<double>&) = 0;
				virtual std::vector<double> getEps() = 0;
				virtual double getMaxEps() = 0;
				virtual std::vector<double> getInnerEps() = 0;
				virtual void printBucketConfig() = 0;
				virtual BucketConfig getConfig() = 0;
				virtual std::vector<common::Triangle<double>> getConvexHull() = 0;
			private:
		};

		class BucketConnected : public Bucket {
		public:
			BucketConnected(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<int, -1, -1>& F);
			std::vector<common::Triangle<double>> getTriangles() override;
			std::vector<common::Edge<double>> getEdges() override;
			std::vector<Eigen::Vector<double, 3>> getVertices() override;
			common::Triangle<double> getTriangle(uint32_t n) override;
			void setEps(std::vector<double>& eps_in) override;
			void setInnerEps(std::vector<double>& inner_eps_in) override;
			std::vector<double> getEps() override;
			std::vector<double> getInnerEps() override;
			void printBucketConfig() override;
			std::vector<common::Triangle<double>> getConvexHull() override;
			BucketConfig getConfig() override;
			double getMaxEps() override;
		private:
			Eigen::Matrix<double, -1, -1> _vertices;
			Eigen::Matrix<int, -1, -1> _faces;
			std::vector<common::Triangle<double>> _hull;
			std::array<double, 8> _eps;
			std::array<double, 8> _inner_eps;
			int _num_tris;
		};

		class BucketSoup : public Bucket {
		public:
			BucketSoup(std::vector<common::Triangle<double>>& tris);
			BucketSoup(const Eigen::Matrix<double, -1, -1>& V, const Eigen::Matrix<int, -1, -1>& F);
			std::vector<common::Triangle<double>> getTriangles() override;
			std::vector<common::Edge<double>> getEdges() override;
			std::vector<Eigen::Vector<double, 3>> getVertices() override;
			common::Triangle<double> getTriangle(uint32_t n) override;
			void setEps(std::vector<double>& eps_in) override;
			void setInnerEps(std::vector<double>& inner_eps_in) override;
			std::vector<double> getEps() override;
			std::vector<double> getInnerEps() override;
			void printBucketConfig() override;
			std::vector<common::Triangle<double>> getConvexHull() override;
			BucketConfig getConfig() override;
			double getMaxEps() override;
		private:
			std::array<common::Triangle<double>, 8> _triangles;
			std::vector<common::Triangle<double>> _hull;
			std::array<double, 8> _eps;
			std::array<double, 8> _inner_eps;
			int _num_tris;
		};
    }
}

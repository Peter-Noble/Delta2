#pragma once

#include "surrogate_tree.h"
#include <Eigen/Dense>

#include "../embree_common/math/bbox.h"

namespace Delta2 {
    class MeshData {
        public:
            MeshData(std::string path, common::Options& opt, bool isSurface = false);
            MeshData(Eigen::MatrixXd V, Eigen::MatrixXi F, common::Options& opt, bool isSurface = false);
            Eigen::Matrix3d getUnitInertiaTensor();
            Eigen::Matrix3d getUnitInertiaTensorInverse();
            double getVolume();
            double getBoundingRadius();

            embree::BBox3fa getAABB() const;
            embree::BBox3fa getAABB(const Eigen::Matrix4d& transform) const;
            embree::BBox3fa getSphereBounds();

            model::SurrogateTree& getSurrogateTree();

            const Eigen::MatrixXd& getVertices();
            const Eigen::MatrixXi& getFaces();

            const std::string serialise();
        private:
            model::SurrogateTree _surrogate;

            Eigen::MatrixXd _vertices;
            Eigen::MatrixXi _faces;

            embree::BBox3fa _AABB;
            double _boundingRadius;
            void initBoudingData();
            void initInertiaDataVolume();
            void initInertiaDataSurface(double thickness);

            double _volume;

            Eigen::Matrix3d _unit_inertia_tensor;
            Eigen::Matrix3d _unit_inertia_tensor_inverse;
    };
}

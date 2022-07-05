#pragma once

#include "igl/opengl/glfw/Viewer.h"

#include "triangle.h"
#include "../model/particle.h"

#include <chrono>
#include <mutex>

namespace Delta2 {
    namespace common {
        class Viewer {
            public:
                Viewer();
                template<typename real>
                void addTriangle(Triangle<real> T) {
                    Eigen::MatrixXd V;
                    Eigen::MatrixXi F;
                    
                    T.toEigen(V, F);

                    concatMesh(_V, _F, V, F, _V, _F);
                }

                void addEigenMesh(const Eigen::Matrix4d& V, const Eigen::MatrixXi& F);
                void addParticle(const Particle& P);
                void addParticleFuture(const Particle& P);
                void addParticleInterval(const Particle& P);
                void addEdge(common::Edge<double> edge);
            
                void show();
                igl::opengl::glfw::Viewer& getViewer();
            private:
                Eigen::MatrixXd _V;
                Eigen::MatrixXi _F;
                igl::opengl::glfw::Viewer _viewer;
                std::vector<common::Edge<double>> _edges;
        };
        class AnimationViewer {
            public:
                AnimationViewer(std::vector<Delta2::Particle>* particles);
                void recordFrame(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& edges, double time);

                AnimationViewer(const AnimationViewer&) = delete;
                void operator=(const AnimationViewer&) = delete;

                void show();
            private:
                std::vector<Delta2::Particle>* _Ps;
                std::vector<std::vector<Eigen::Matrix4d>> _Ts;
                std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> _edges;
                std::vector<double> _times;
                std::vector<std::vector<Eigen::Vector3d>> _colours;
                igl::opengl::glfw::Viewer _viewer;
                int _frame;
                std::chrono::time_point<std::chrono::system_clock> _start_time;
                std::mutex _lock;
        };
    }
}

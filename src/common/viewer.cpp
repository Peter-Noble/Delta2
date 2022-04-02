#include "viewer.h"
#include "basic_utils.h"
#include <chrono>

using namespace Delta2::common;

Viewer::Viewer() {
    _viewer.data().face_based = true;
}

void Viewer::show() {
    Eigen::Matrix4d T = Delta2::common::transformationMatrix(Eigen::Vector3d({-igl::PI/2.0, 0.0, 0.0}), Eigen::Vector3d({0.0, 0.0, 0.0}));
    Eigen::MatrixXd V = transform(_V, T);

    _viewer.data().set_mesh(V, _F);

    _viewer.launch(true, false, "Delta2", 0, 0);
}

void Viewer::addEigenMesh(const Eigen::Matrix4d& V, const Eigen::MatrixXi& F) {
    concatMesh(_V, _F, V, F, _V, _F);
}

void Viewer::addParticle(const Particle& P) {
    concatMesh(_V, _F, P.getTransformedVertices(), P.getFaces(), _V, _F);
}

void Viewer::addParticleInterval(const Particle& P) {
    concatMesh(_V, _F, P.getTransformedVertices(), P.getFaces(), _V, _F);
    concatMesh(_V, _F, P.getTransformedVerticesFuture(), P.getFaces(), _V, _F);
}

igl::opengl::glfw::Viewer& Viewer::getViewer() {
    return _viewer;
}

AnimationViewer::AnimationViewer(std::vector<Delta2::Particle>* particles) : _Ps(particles) {
}

void AnimationViewer::recordFrame(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& edges, double time) {
    std::lock_guard<std::mutex> guard(_lock);
    std::vector<Eigen::Matrix4d> frame_transforms;
    std::vector<Eigen::Vector3d> colours;
    // _Ts.emplace_back();
    for (Delta2::Particle& p : *_Ps) {
        // _Ts[_Ts.size() - 1].push_back(p.current_state.getTransformation());
        frame_transforms.push_back(p.current_state.getTransformation());
        if (p.is_static) {
            colours.push_back({0.3, 0.3, 0.4});
        }
        else if (p.getSleeping()) {
            colours.push_back({0.5, 0.4, 0.1});
        }
        else {
            colours.push_back({1.0, 0.7, 0.3});
        }
    }
    _Ts.push_back(frame_transforms);
    _edges.push_back(edges);
    _times.push_back(time);
    _colours.push_back(colours);
}

void AnimationViewer::show() {
    _viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool {
        std::lock_guard<std::mutex> guard(_lock);
        if (_viewer.core().is_animating && _Ts.size() > 0) {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::MatrixXd C;

            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = current_time - _start_time;

            double seconds = std::chrono::duration<double>(elapsed).count();

            _frame = 0;
            // printf("times: %f, seconds: %f\n", _times[_frame], seconds);
            while (_times[_frame] < seconds) {
                _frame++;
                if (_frame >= _times.size()) {
                    _start_time = current_time;
                    seconds = 0.0;
                    _frame = 0;
                }
            }

            for (int p_i = 0; p_i < (*_Ps).size(); p_i++) {
                Particle& p = (*_Ps)[p_i];
                Eigen::MatrixXi Fp = p.mesh->getFaces();
                Eigen::MatrixXd Vp = transform(p.mesh->getVertices(), _Ts[_frame][p_i]);
                concatMesh(V, F, Vp, Fp, V, F);
                Eigen::MatrixXd Cp;
                Cp.resize(Fp.rows(), 3);
                for (int i = 0; i < Fp.rows(); i++) {
                    for (int j = 0; j < 3; j++) {
                        Cp(i,j) = _colours[_frame][p_i][j];
                    }
                }
                concatColors(C, Cp, C);
            }

            Eigen::Matrix4d T = transformationMatrix<double>({-igl::PI/2.0, 0, 0}, {0, 0, 0});
            V = transform(V, T);
            _viewer.data().set_mesh(V, F);
            _viewer.data().set_colors(C);

            int num_edges = _edges[_frame].size();
            Eigen::MatrixXd a_pts;
            a_pts.resize(num_edges, 3);
            Eigen::MatrixXd b_pts;
            b_pts.resize(num_edges, 3);
            for (int e_i = 0; e_i < num_edges; e_i++) {
                a_pts.row(e_i) = transform(_edges[_frame][e_i].first, T);
                b_pts.row(e_i) = transform(_edges[_frame][e_i].second, T);
            }
            _viewer.data().clear_points();
            _viewer.data().clear_edges();
            _viewer.data().add_points(a_pts, Eigen::RowVector3d(1.0, 0.0, 0.0));
            _viewer.data().add_points(b_pts, Eigen::RowVector3d(0.0, 0.0, 1.0));
            _viewer.data().add_edges(a_pts, b_pts, Eigen::RowVector3d(1.0, 1.0, 1.0));
        }
        return false;
    };
    _viewer.core().is_animating = true;
    _frame = 0;
    _start_time = std::chrono::high_resolution_clock::now();
    _viewer.data().set_face_based(true);

    _viewer.launch(true, false, "Delta2", 0, 0);
}
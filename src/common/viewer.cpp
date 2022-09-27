#include "viewer.h"
#include "basic_utils.h"
#include <chrono>
#include <tuple>

using namespace Delta2::common;

Viewer::Viewer() {
    _viewer.data().face_based = true;
}

void Viewer::show() {
    printf("Viewer show\n");
    Eigen::Matrix4d T = Delta2::common::transformationMatrix(Eigen::Vector3d({-igl::PI/2.0, 0.0, 0.0}), Eigen::Vector3d({0.0, 0.0, 0.0}));
    Eigen::MatrixXd V = transform(_V, T);

    _viewer.data().set_mesh(V, _F);

    int num_edges = _edges.size();
    Eigen::MatrixXd a_pts;
    a_pts.resize(num_edges, 3);
    Eigen::MatrixXd b_pts;
    b_pts.resize(num_edges, 3);
    for (int e_i = 0; e_i < num_edges; e_i++) {
        a_pts.row(e_i) = transform(_edges[e_i].A, T);
        b_pts.row(e_i) = transform(_edges[e_i].B, T);
    }
    _viewer.data().add_points(a_pts, Eigen::RowVector3d(1.0, 0.0, 0.0));
    _viewer.data().add_points(b_pts, Eigen::RowVector3d(0.0, 0.0, 1.0));
    _viewer.data().add_edges(a_pts, b_pts, Eigen::RowVector3d(1.0, 1.0, 1.0));

    _viewer.launch(true, false, "Delta2", 0, 0);
}

void Viewer::addEigenMesh(const Eigen::Matrix4d& V, const Eigen::MatrixXi& F) {
    concatMesh(_V, _F, V, F, _V, _F);
}

void Viewer::addParticle(const Particle& P) {
    concatMesh(_V, _F, P.getTransformedVertices(), P.getFaces(), _V, _F);
}

void Viewer::addParticleFuture(const Particle& P) {
    concatMesh(_V, _F, P.getTransformedVerticesFuture(), P.getFaces(), _V, _F);
}

void Viewer::addParticleInterval(const Particle& P) {
    concatMesh(_V, _F, P.getTransformedVertices(), P.getFaces(), _V, _F);
    concatMesh(_V, _F, P.getTransformedVerticesFuture(), P.getFaces(), _V, _F);
}

void Viewer::addEdge(common::Edge<double> edge) {
    _edges.push_back(edge);
}

igl::opengl::glfw::Viewer& Viewer::getViewer() {
    return _viewer;
}

AnimationViewer::AnimationViewer(std::vector<Delta2::Particle>* particles) : _Ps(particles) {
}

void AnimationViewer::recordFrame(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& edges) {
    std::lock_guard<std::mutex> guard(_lock);
    std::vector<Eigen::Matrix4d> frame_transforms;
    std::vector<Eigen::Vector3d> colours;
    std::vector<int> groups;
    std::vector<ParticleState> states;
    // _Ts.emplace_back();

    double time = std::numeric_limits<double>::infinity();
    for (Delta2::Particle& p : *_Ps) {
        if (!p.is_static) {
            time = std::min(time, p.current_state.getTime());
        }
    }
    if (std::isinf(time)) {
        time = 0.0;
    }

    for (Delta2::Particle& p : *_Ps) {
        // _Ts[_Ts.size() - 1].push_back(p.current_state.getTransformation());
        if (p.is_static) {
            frame_transforms.push_back(p.current_state.getTransformation());
        }
        else {
            frame_transforms.push_back(p.current_state.interpolate(p.last_state, time).getTransformation());
        }
        if (p.is_static) {
            colours.push_back({0.3, 0.3, 0.4});
            states.push_back(ParticleState::Static);
        }
        else if (p.getSleeping()) {
            colours.push_back({0.5, 0.4, 0.1});
            states.push_back(ParticleState::Sleeping);
        }
        else {
            colours.push_back({1.0, 0.7, 0.3});
            states.push_back(ParticleState::Active);
        }
        groups.push_back(p.cluster_id);
    }
    _edges.push_back(edges);
    _times.push_back(time);
    _colours.push_back(colours);
    _states.push_back(states);
    _groups.push_back(groups);
    _Ts.push_back(frame_transforms);
    _contacts.push_back({});
}

void AnimationViewer::recordFrame(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>>& contacts) {
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> edges;
    recordFrame(edges);
    _contacts[_contacts.size() - 1] = contacts;
}

void AnimationViewer::show() {
    printf("Animation show\n");
    std::vector<Eigen::Vector3d> group_colours {
        {230.0 / 255.0, 25.0 / 255.0,   75.0 / 255.0},
        {60.0 / 255.0,  180.0 / 255.0,  75.0 / 255.0},
        {255.0 / 255.0, 255.0 / 255.0,  25.0 / 255.0},
        {0.0 / 255.0,   130.0 / 255.0,  200.0 / 255.0},
        {245.0 / 255.0, 130.0 / 255.0,  48.0 / 255.0},
        {145.0 / 255.0, 30.0 / 255.0,   180.0 / 255.0},
        {70.0 / 255.0,  240.0 / 255.0,  240.0 / 255.0},
        {240.0 / 255.0, 50.0 / 255.0,   230.0 / 255.0},
        {210.0 / 255.0, 245.0 / 255.0,  60.0 / 255.0},
        {250.0 / 255.0, 190.0 / 255.0,  212.0 / 255.0},
        {0.0 / 255.0,   128.0 / 255.0,  128.0 / 255.0},
        {220.0 / 255.0, 190.0 / 255.0,  255.0 / 255.0},
        {170.0 / 255.0, 110.0 / 255.0,  40.0 / 255.0},
        {255.0 / 255.0, 250.0 / 255.0,  200.0 / 255.0},
        {128.0 / 255.0, 0.0 / 255.0,    0.0 / 255.0},
        {170.0 / 255.0, 255.0 / 255.0,  195.0 / 255.0},
        {128.0 / 255.0, 128.0 / 255.0,  0.0 / 255.0},
        {255.0 / 255.0, 215.0 / 255.0,  180.0 / 255.0},
        {0.0 / 255.0,   0.0 / 255.0,    128.0 / 255.0}
    };

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
                Eigen::Vector3d use_colour;
                assert(_frame < _states.size());
                assert(p_i < _states[_frame].size());
                switch (_states[_frame][p_i]) {
                    case ParticleState::Active:
                        use_colour = group_colours[_groups[_frame][p_i] % group_colours.size()];
                        break;
                    case ParticleState::Sleeping:
                        use_colour = {0.5, 0.4, 0.1};
                        break;
                    case ParticleState::Static:
                        use_colour = {0.3, 0.3, 0.4};
                        break;
                }
                for (int i = 0; i < Fp.rows(); i++) {
                    for (int j = 0; j < 3; j++) {
                        Cp(i,j) = use_colour[j];
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

            Eigen::MatrixXd hit_normal_edges;
            hit_normal_edges.resize(_contacts[_frame].size(), 3);
            Eigen::MatrixXd hit_friction_edges;
            hit_friction_edges.resize(_contacts[_frame].size(), 3);
            Eigen::MatrixXd hit_pts;
            hit_pts.resize(_contacts[_frame].size(), 3);
            for (int c_i = 0; c_i < _contacts[_frame].size(); c_i++) {
                hit_pts.row(c_i) = transform(std::get<0>(_contacts[_frame][c_i]), T);
                hit_normal_edges.row(c_i) = transform(Eigen::Vector3d(std::get<0>(_contacts[_frame][c_i]) + std::get<1>(_contacts[_frame][c_i])), T);
                hit_friction_edges.row(c_i) = transform(Eigen::Vector3d(std::get<0>(_contacts[_frame][c_i]) + std::get<2>(_contacts[_frame][c_i])), T);
            }

            _viewer.data().clear_points();
            _viewer.data().clear_edges();
            _viewer.data().add_points(a_pts, Eigen::RowVector3d(1.0, 0.0, 0.0));
            _viewer.data().add_points(b_pts, Eigen::RowVector3d(0.0, 0.0, 1.0));
            _viewer.data().add_edges(a_pts, b_pts, Eigen::RowVector3d(1.0, 1.0, 1.0));

            _viewer.data().add_points(hit_pts, Eigen::RowVector3d(1.0, 1.0, 0.0));
            _viewer.data().add_edges(hit_pts, hit_friction_edges, Eigen::RowVector3d(0.0, 1.0, 0.0));
            _viewer.data().add_edges(hit_pts, hit_normal_edges, Eigen::RowVector3d(0.0, 1.0, 1.0));

        }
        return false;
    };
    _viewer.core().is_animating = true;
    _frame = 0;
    _start_time = std::chrono::high_resolution_clock::now();
    _viewer.data().set_face_based(true);

    _viewer.launch(true, false, "Delta2", 0, 0);
}
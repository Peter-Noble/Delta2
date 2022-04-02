#pragma once

#include "mesh.h"

#include <memory>

namespace Delta2 {
    class State {
    public:
        State();
        double getTime() const;
        void setTime(double t);
        const Eigen::Quaterniond& getRotation() const;
        void setRotation(const Eigen::Quaterniond& r);
        const Eigen::Vector3d& getAngularMomentum() const;
        void setAngular(const Eigen::Vector3d& momentum);
        const Eigen::Vector3d& getTranslation() const;
        void setTranslation(Eigen::Vector3d t);
        const Eigen::Vector3d& getVelocity() const;
        void setVelocity(Eigen::Vector3d v);
        Eigen::Matrix4d getTransformation() const;

        State extrapolate(double t, Eigen::Matrix3d inv_inertia);
        Eigen::Vector3d pointVelocity(const Eigen::Vector3d& pt, const State& future) const;
        bool isValid() const;
        State interpolate(State last, double time) const;
    private:
        double _time;

        Eigen::Quaterniond _rotation;
        Eigen::Vector3d _angular_momentum;
        Eigen::Vector3d _translation;
        Eigen::Vector3d _velocity;
    };

    class Particle {
    public:
        Particle(std::shared_ptr<MeshData> mesh, double density, double friction_coeff, double eps);
        Particle(std::shared_ptr<MeshData> mesh, double density, double friction_coeff, double eps, State current, State future);
        void operator=(const Particle&) = delete;

        Particle copy();

        double getMass() const;

        Eigen::MatrixXd getTransformedVertices() const;
        Eigen::MatrixXd getTransformedVerticesFuture() const;
        Eigen::MatrixXi getFaces() const;

        void setSleeping(bool s);
        bool getSleeping() const;

        Eigen::Matrix3d getInverseInertiaMatrix() const;

        void projectFutureState(double t);

        Eigen::Vector3d futurePointVelocity(const Eigen::Vector3d& pt) const;

        std::shared_ptr<MeshData> mesh;

        void rollBackState(double time);

        State last_state;
        State current_state;
        State future_state;

        double maxDifferenceToCurrent(const State& S);

        bool is_static;
        double geo_eps;
        double friction_coeff;
        double last_time_step_size;
        int id;
        double sleep_candidate_time;
    private:
        void assignID();

        double _density;
        double _mass;
        bool _is_sleeping;

        bool _got_inertia_body;
        Eigen::Matrix3d _inertia_body;
        Eigen::Matrix3d _inertia_body_inverse;
    };
}

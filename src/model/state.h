#pragma once

#include <Eigen/Dense>

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
            Eigen::Vector3d pointVelocity(const Eigen::Vector3d& pt, const Eigen::Matrix3d& inv_inertia) const;
            bool isValid() const;
            State interpolate(State last, double time) const;
            bool isStationary() const;

            void applyDelta(double t, Eigen::Vector3d F, Eigen::Vector3d T, double mass, const Eigen::Matrix3d& inverse_inertia);

            std::string serialise();
        private:
            double _time;

            Eigen::Quaterniond _rotation;
            Eigen::Vector3d _angular_momentum;
            Eigen::Vector3d _translation;
            Eigen::Vector3d _velocity;
    };
}

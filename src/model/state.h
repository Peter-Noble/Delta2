#pragma once

#include <Eigen/Dense>
#include "../common/utils.h"

namespace Delta2 {
    class State {
        public:
            State();
            inline double getTime() const {return _time;};
            inline void setTime(double t) {_time = t;};
            inline const Eigen::Quaterniond& getRotation() const {return _rotation;};
            inline void setRotation(const Eigen::Quaterniond& r) {_rotation = r;};
            inline const Eigen::Vector3d& getAngularMomentum() const {return _angular_momentum;};
            inline void setAngular(const Eigen::Vector3d& momentum) {_angular_momentum = momentum;};
            inline const Eigen::Vector3d& getTranslation() const {return _translation;};
            inline void setTranslation(Eigen::Vector3d t) {_translation = t;};
            inline const Eigen::Vector3d& getVelocity() const {return _velocity;};
            inline void setVelocity(Eigen::Vector3d v) {_velocity = v;};
            Eigen::Matrix4d getTransformation() const {return common::transformationMatrix(_rotation, _translation);};

            State extrapolate(double t, Eigen::Matrix3d inv_inertia);
            Eigen::Vector3d pointVelocity(const Eigen::Vector3d& pt, const State& future) const {
                Eigen::Vector4d new_pt = future.getTransformation() * (getTransformation().inverse() * pt.homogeneous());
                return (new_pt - pt.homogeneous()).head<3>() / (future.getTime() - getTime());
            };
            Eigen::Vector3d pointVelocity(const Eigen::Vector3d& pt, const Eigen::Matrix3d& inv_inertia) const {
                Eigen::Vector3d w = inv_inertia * getAngularMomentum();
                return getVelocity() + w.cross(pt - getTranslation());
            };
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

#include "state.h"
#include "../common/utils.h"
#include "../globals.h"

using namespace Delta2;
using namespace Eigen;

State::State() {
    _time = 0;
    _rotation = Eigen::Quaterniond::Identity();
    _angular_momentum = { 0, 0, 0 };
    _translation = { 0, 0, 0 };
    _velocity = { 0, 0, 0 };
}

double State::getTime() const {
    return _time;
}

void State::setTime(double t) {
    _time = t;
}

const Eigen::Quaterniond& State::getRotation() const {
    return _rotation;
}

void State::setRotation(const Eigen::Quaterniond& r) {
    _rotation = r;
}

const Eigen::Vector3d& State::getAngularMomentum() const {
    return _angular_momentum;
}

void State::setAngular(const Eigen::Vector3d& momentum) {
    _angular_momentum = momentum;
}

const Eigen::Vector3d& State::getTranslation() const {
    return _translation;
}

void State::setTranslation(Eigen::Vector3d t) {
    _translation = t;
}

const Eigen::Vector3d& State::getVelocity() const {
    return _velocity;
}

void State::setVelocity(Eigen::Vector3d v) {
    _velocity = v;
}

State State::interpolate(State last, double time) const {
    #ifndef NDEBUG
    if (time > _time + 1e-4 || time + 1e-4 < last.getTime() || last.getTime() > _time + 1e-4) {
        throw std::runtime_error("Invalid args for interpolation");
    }
    #endif
    if (last.getTime() == _time) {
        return *this;
    }
    double t = (time - last.getTime()) / (_time - last.getTime());
    State result;
    result.setTime(time);
    result.setTranslation(t * _translation + (1-t) * last.getTranslation());
    result.setVelocity(t * _velocity + (1-t) * last.getVelocity());
    result.setAngular(t * _angular_momentum + (1-t) * last.getAngularMomentum());
    result.setRotation(last.getRotation().slerp(t, _rotation));
    return result;
}

// pt is in global space
// the returned value is in relative space to the current orientation in this state
Eigen::Vector3d State::pointVelocity(const Eigen::Vector3d& pt, const State& future) const {
    // TODO velocity at a point section - https://gafferongames.com/post/physics_in_3d/
    Eigen::Vector4d new_pt = future.getTransformation() * (getTransformation().inverse() * pt.homogeneous());
    return (new_pt - pt.homogeneous()).head<3>() / (future.getTime() - getTime());
}

// pt is in global space
Eigen::Vector3d State::pointVelocity(const Eigen::Vector3d& pt, const Eigen::Matrix3d& inv_inertia) const {
    Eigen::Vector3d w = inv_inertia * getAngularMomentum();
    // globals::logger.printf(2, "w: %f, %f, %f\n", w.x(), w.y(), w.z());
    return getVelocity() + w.cross(pt - getTranslation());
}

Eigen::Matrix4d State::getTransformation() const {
    return common::transformationMatrix(_rotation, _translation);
}

State State::extrapolate(double t, Eigen::Matrix3d inv_inertia) {
    State result;
    result.setTime(_time + t);
    result.setVelocity(getVelocity());
    result.setTranslation(getTranslation() + t * getVelocity());
    Eigen::Matrix3d R = _rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * inv_inertia * R.transpose();
    Eigen::Vector3d omega = Iinv * _angular_momentum;
    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);
    result.setRotation(_rotation * rot_delta);
    result.setAngular(getAngularMomentum());
    return result;
}

bool State::isValid() const {
    // if (_translation.hasNaN()) {
    //     return false;
    // }
    // if (_velocity.hasNaN()) {
    //     return false;
    // }
    // if (_rotation.coeffs().hasNaN()) {
    //     return false;
    // }
    // if (_angular_momentum.hasNaN()) {
    //     return false;
    // }
    return true;
}

bool State::isStationary() const {
    return getVelocity().norm() < 1e-3 && getAngularMomentum().norm() < 1e-2;
}

void State::applyDelta(double t, Eigen::Vector3d force, Eigen::Vector3d torque, double mass, const Eigen::Matrix3d& inverse_inertia) {
    setTime(getTime() + t);
    setVelocity(getVelocity() + t * force / mass);
    setTranslation(getTranslation() + t * getVelocity());

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = getAngularMomentum() + t * torque;

    Eigen::Quaterniond rotation = getRotation();

    setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * inverse_inertia * R.transpose();
    Eigen::Vector3d omega = Iinv * angular_momentum;
    
    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;

    setRotation(rotation);
}

std::string State::serialise() {
    std::stringstream result;
    result << _time << "\n";
    // Eigen::Vector3d r = _rotation.toRotationMatrix().eulerAngles(0, 1, 2);
    // result << r.x() << " " << r.y() << " " << r.z() << " " << "\n";
    result << _rotation.w() << " " << _rotation.x() << " " << _rotation.y() << " " << _rotation.z() << " " << "\n";
    // result << _angular_momentum.x() << " " << _angular_momentum.y() << " " << _angular_momentum.z() << "\n";
    result << _translation.x() << " " << _translation.y() << " " << _translation.z() << "\n";
    // result << _velocity.x() << " " << _velocity.y() << " " << _velocity.z() << "\n";
    return result.str();
}
#include "particle.h"
#include "../common/basic_utils.h"

using namespace Delta2;

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
    if (time > _time || time < last.getTime() || last.getTime() > _time) {
        throw std::runtime_error("Invalid args for interpolation");
    }
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

// Particle

Particle::Particle(std::shared_ptr<MeshData> m, double density, double friction_co, double eps) {
    mesh = m;
    _inertia_body = mesh->getUnitInertiaTensor() * density;
    _inertia_body_inverse = mesh->getUnitInertiaTensorInverse() * 1.0/density;
    _density = density;
    _mass = _density * mesh->getVolume();
    friction_coeff = friction_co;
    is_static = false;
    geo_eps = eps;
    last_time_step_size = 1.0;
    _is_sleeping = false;
    sleep_candidate_time = 0.0;
    assignID();
}

Particle::Particle(std::shared_ptr<MeshData> m, double density, double friction_co, double eps, State current, State future) {
    mesh = m;
    _inertia_body = mesh->getUnitInertiaTensor() * density;
    _inertia_body_inverse = mesh->getUnitInertiaTensorInverse() * 1.0/density;
    _density = density;
    _mass = _density * mesh->getVolume();
    friction_coeff = friction_co;
    current_state = current;
    future_state = future;
    is_static = false;
    geo_eps = eps;
    last_time_step_size = 1.0;
    _is_sleeping = false;
    assignID();
}

void Particle::assignID() {
    static uint32_t ID_counter = 0;
    id = ID_counter++;
}

double Particle::getMass() const {
    return _mass;
}

Eigen::MatrixXd Particle::getTransformedVertices() const {
    return common::transform(mesh->getVertices(), current_state.getTransformation());
}

Eigen::MatrixXd Particle::getTransformedVerticesFuture() const {
    return common::transform(mesh->getVertices(), future_state.getTransformation());
}

Eigen::MatrixXi Particle::getFaces() const {
    return mesh->getFaces();
}

void Particle::projectFutureState(double t) {
    future_state = current_state.extrapolate(t, _inertia_body_inverse);
}

Eigen::Matrix3d Particle::getInverseInertiaMatrix() const {
    return _inertia_body_inverse;
}

Eigen::Vector3d Particle::futurePointVelocity(const Eigen::Vector3d& pt) const {
    Eigen::Matrix3d R = future_state.getRotation().toRotationMatrix();
    Eigen::Matrix3d Iinv = R * _inertia_body_inverse * R.transpose();
    Eigen::Vector3d omega = Iinv * future_state.getAngularMomentum();

    return future_state.getVelocity() + omega.cross(pt - future_state.getTranslation());
}

void Particle::rollBackState(double time) {
    current_state = current_state.interpolate(last_state, time);
}

void Particle::setSleeping(bool s) {
    _is_sleeping = s;
}

bool Particle::getSleeping() const {
    return _is_sleeping || is_static;
}

double Particle::maxDifferenceToCurrent(const State& S) {
    double ang = S.getRotation().angularDistance(current_state.getRotation());
    double rot_dist = mesh->getBoundingRadius() * std::tan(ang);
    return (rot_dist + (current_state.getTranslation() - S.getTranslation()).norm()) / (S.getTime() - current_state.getTime());
}

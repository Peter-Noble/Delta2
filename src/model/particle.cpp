#include "particle.h"
#include "../common/basic_utils.h"

using namespace Delta2;
using namespace Eigen;

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
    sleep_not_moved_steps = 0;
    restitution = 0.5;
    cluster_id = -1;
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
    sleep_not_moved_steps = 0;
    _is_sleeping = false;
    cluster_id = -1;
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

Eigen::Vector3d Particle::pointVelocity(const Eigen::Vector3d& pt) const {
    Eigen::Matrix3d R = current_state.getRotation().toRotationMatrix();
    Eigen::Matrix3d Iinv = R * _inertia_body_inverse * R.transpose();
    Eigen::Vector3d omega = Iinv * current_state.getAngularMomentum();

    return current_state.getVelocity() + omega.cross(pt - current_state.getTranslation());
}

Eigen::Vector3d Particle::pointVelocity(const Eigen::Vector3d& pt, const State& S) const {
    Eigen::Matrix3d R = S.getRotation().toRotationMatrix();
    Eigen::Matrix3d Iinv = R * _inertia_body_inverse * R.transpose();
    Eigen::Vector3d omega = Iinv * S.getAngularMomentum();

    return S.getVelocity() + omega.cross(pt - S.getTranslation());
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

State Particle::updateState(double t, const Vector3d& force, const Vector3d& torque) {
    State updated = current_state;

    updated.setTime(current_state.getTime() + t);
    Vector3d vel_change = t * force / getMass();
    updated.setVelocity(current_state.getVelocity() + t * force / getMass());
    // updated.setTranslation(getTranslation() + t * updated.getVelocity() + force_offset);
    Vector3d pos_change_from_velocity = t * current_state.getVelocity();
    Vector3d pos_change_from_accel = 0.5 * force / getMass() * t * t;
    updated.setTranslation(current_state.getTranslation() + t * current_state.getVelocity() + 0.5 * force / getMass() * t * t);

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = current_state.getAngularMomentum() + t * torque;

    Eigen::Quaterniond rotation = current_state.getRotation();

    updated.setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * getInverseInertiaMatrix() * R.transpose();
    Eigen::Vector3d omega = Iinv * angular_momentum;

    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;
    updated.setRotation(rotation);

    return updated;
}

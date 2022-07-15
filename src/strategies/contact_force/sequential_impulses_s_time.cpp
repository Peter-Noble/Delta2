#include "sequential_impulses_s_time.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"

#include "../../common/viewer.h"

using namespace Delta2;
using namespace strategy;
using namespace Eigen;

struct LocalContactS {
    Eigen::Vector3d local_A;
    Eigen::Vector3d local_B;
    Eigen::Vector3d global_centre;
    Eigen::Vector3d global_normal;
    Eigen::Vector3d global_tangent;
    double mass_eff_normal;
    double mass_eff_tangent;
    double s;
    State S_A;
    State S_B;

    double next_primary_score;

    double force_last;
    double force_move_offset;
};

double averageForce(double f0, double f1, double d0, double d1, double outer) {
    assert(f0 >= 0.0);
    assert(f1 >= 0.0);
    assert(outer > 0.0);
    return (f0 + f1) / 2.0;
    // if (d0 > outer && d1 > outer) {
    //     return 0.0;
    // }
    // else if (d0 > outer) {
    //     double t = (outer - d0) / (d1 - d0); // the point at which d == outer
    //     return (f0 + f1) * (1.0 - t) / 2.0;
    // }
    // else if (d1 > outer) {
    //     double t = (outer - d0) / (d1 - d0); // the point at which d == outer
    //     return (f0 + f1) * t / 2.0;
    // }
    // else {
    //     return (f0 + f1) / 2.0;
    // }
}

double forceDamp(double force, double v0, double v1, double damp) {
    double u = std::clamp(v1 / (v1 - v0), 0.0, 1.0);
    return force * (1.0 - u) + force * u * damp;
}

State updateState(State& current, double t, const Vector3d& force, const Vector3d& force_offset, const Vector3d& torque, const Vector3d& torque_offset, Particle* p) {
    State updated = current;

    updated.setTime(current.getTime() + t);
    Vector3d vel_change = t * force / p->getMass();
    updated.setVelocity(current.getVelocity() + t * force / p->getMass());
    // updated.setTranslation(current.getTranslation() + t * updated.getVelocity() + force_offset);
    Vector3d pos_change_from_velocity = t * current.getVelocity();
    Vector3d pos_change_from_accel = 0.5 * force / p->getMass() * t * t;
    updated.setTranslation(current.getTranslation() + t * current.getVelocity() + 0.5 * force / p->getMass() * t * t);

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = current.getAngularMomentum() + t * torque;

    Eigen::Quaterniond rotation = current.getRotation();

    updated.setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * p->getInverseInertiaMatrix() * R.transpose();
    Eigen::Vector3d omega = Iinv * (angular_momentum + t * t * torque_offset);

    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;
    updated.setRotation(rotation);

    return updated;
}

State updateStateWithOffset(State& current, double t, const Vector3d& force, const Vector3d& force_offset, const Vector3d& torque, const Vector3d& torque_offset, Particle* p) {
    State updated = current;

    updated.setTime(current.getTime() + t);
    updated.setVelocity(current.getVelocity() + t * force / p->getMass());
    updated.setTranslation(current.getTranslation() + t * t * force_offset / p->getMass());

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = current.getAngularMomentum() + t * torque;

    Eigen::Quaterniond rotation = current.getRotation();

    updated.setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * p->getInverseInertiaMatrix() * R.transpose();
    Eigen::Vector3d omega = Iinv * (angular_momentum + t * t * torque_offset);

    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;
    updated.setRotation(rotation);

    return updated;
}

SequentialImpulsesSTime::SequentialImpulsesSTime(FrictionStrategy& friction,
                                     common::Options& opt) :
                                     ContactForceStrategy(friction, opt) {

}

bool SequentialImpulsesSTime::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) {
    std::vector<LocalContactS> contacts;

    std::sort(hits.begin(), hits.end());

    bool off_normal = false;

    // common::Viewer view;
    // view.addParticleInterval(cluster.particles[0]);
    // view.addParticleInterval(cluster.particles[1]);
    // for (int c = 0; c < hits.size(); c++) {
    //     view.addEdge(common::Edge(hits[c].A, hits[c].B));
    // }
    // view.show();

    for (int c = 0; c < hits.size(); c++) {
        LocalContactS lc;

        Eigen::Matrix4d A_T_i = hits[c].p_a->future_state.getTransformation().inverse();
        lc.local_A = common::transform(hits[c].A, A_T_i);
        Eigen::Matrix4d B_T_i = hits[c].p_b->future_state.getTransformation().inverse();
        lc.local_B = common::transform(hits[c].B, B_T_i);

        lc.global_centre = (hits[c].A + hits[c].B) / 2.0;
        lc.global_normal = (hits[c].B - hits[c].A) / 2.0;
        lc.global_tangent = {0, 0, 0}; // TODO

        if (std::fabs(lc.global_normal.x()) > 1e-8) {
            off_normal = true;
        }
        printf("Normal: %f, %f, %f\n", lc.global_normal.x(), lc.global_normal.y(), lc.global_normal.z());
        printf("A: %f, %f, %f\n", hits[c].A.x(), hits[c].A.y(), hits[c].A.z());
        printf("B: %f, %f, %f\n", hits[c].B.x(), hits[c].B.y(), hits[c].B.z());

        if (lc.global_normal.norm() < 1e-8) {
            common::Viewer view;
            view.addParticleInterval(*hits[c].p_a);
            view.addParticleInterval(*hits[c].p_b);
            view.addEdge(common::Edge(hits[c].A, hits[c].B));
            view.show();

            return false;
            assert(false);
        }

        Eigen::Vector3d r_A = hits[c].A - hits[c].p_a->future_state.getTranslation();
        Eigen::Vector3d r_B = hits[c].B - hits[c].p_b->future_state.getTranslation();

        Eigen::Vector3d rn_A = r_A.cross(lc.global_normal);
        Eigen::Vector3d rn_B = r_B.cross(lc.global_normal);

        Eigen::Vector3d rt_A = r_A.cross(lc.global_tangent);
        Eigen::Vector3d rt_B = r_B.cross(lc.global_tangent);

        double im_A = 1.0 / hits[c].p_a->getMass();
        double im_B = 1.0 / hits[c].p_b->getMass();
        Eigen::Matrix3d ii_A = hits[c].p_a->getInverseInertiaMatrix();
        Eigen::Matrix3d ii_B = hits[c].p_b->getInverseInertiaMatrix();

        double k_normal = 0.0;
        double k_tangent = 0.0;
        if (!hits[c].p_a->is_static) {
            k_normal += im_A + rn_A.transpose() * ii_A * rn_A;
            k_tangent += im_A + rt_A.transpose() * ii_A * rt_A;
            // k_normal += im_A; // TODO
            // k_tangent += im_A;  // TODO

        }
        if (!hits[c].p_b->is_static) {
            k_normal += im_B + rn_B.transpose() * ii_B * rn_B;
            k_tangent += im_B + rt_B.transpose() * ii_B * rt_B;
            // k_normal += im_B;  // TODO
            // k_tangent += im_B;  // TODO
        }

        lc.mass_eff_normal = k_normal < 1e-6 ? 0.0 : 1.0 / k_normal;
        lc.mass_eff_tangent = k_tangent < 1e-6 ? 0.0 : 1.0 / k_tangent;

        lc.s = 1.0;
        lc.S_A = hits[c].p_a->future_state;
        lc.S_B = hits[c].p_b->future_state;

        lc.force_last = 0.0;
        lc.force_move_offset = 0.0;

        contacts.push_back(lc);
    }

    // if (off_normal) {
    //     common::Viewer view;
    //     for (Particle* p : cluster.particles) {
    //         view.addParticleFuture(*p);
    //     }
    //     for (int c = 0; c < hits.size(); c++) {
    //         view.addEdge(common::Edge(hits[c].A, hits[c].B));
    //     }
    //     view.show();
    // }

    std::vector<Eigen::Vector3d> impulses;
    std::vector<Eigen::Vector3d> rotational_impulses;
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<Eigen::Vector3d> forces_last;
    std::vector<Eigen::Vector3d> torques_last;
    std::vector<Eigen::Vector3d> force_offsets;
    std::vector<Eigen::Vector3d> torque_offsets;
    std::vector<double> min_s;
    min_s.resize(cluster.particles.size());

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        forces.push_back(external_force(cluster.particles[p_i]));
        torques.push_back({0, 0, 0});
        forces_last.push_back(external_force(cluster.particles[p_i]));
        torques_last.push_back({0, 0, 0});
        force_offsets.push_back({0, 0, 0});
        torque_offsets.push_back({0, 0, 0});
    }

    double last_ds;

    double damp_start = 0.75;

    std::vector<State> SStates;
    std::vector<State> FStates;
    for (int p = 0; p < cluster.particles.size(); p++) {
        SStates.push_back(cluster.particles[p].current_state);
        FStates.push_back(cluster.particles[p].future_state);
        min_s[p] = 1.0;
    }

    int pre_S_iterations = 200;
    int s_calc_iterations = 100;
    int new_S_iterations = 100;
    int contact_point_s_calc_iterations = 100;
    int post_S_iterations = 1000;
    int post_S_distance_iterations = 100;

    double s_diff = 0.2;
    double s_diff_adjust = 0.91;

    double force_damp = 1.0;
    double force_damp_adjust = 0.95;
    
    double separation_damp = 0.8;  // The proportion of energy used for a contact when it's separating

    double s_damp = 0.1;
    int primary = 0;
    double s = 0.0;
    if (hits.size()) {
        for (int it = 0; it < pre_S_iterations; it++) {
            // printf("Iteration: %i\n", it);
            double s_last = s;

            // Select a new timestep
            int primary_a_id, primary_b_id;
            {
                uint32_t a_id = cluster.particles.getLocalID(hits[primary].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[primary].p_b->id);
                primary_a_id = a_id;
                primary_b_id = b_id;

                if (!hits[primary].p_a->is_static)
                {
                    Eigen::Vector3d force = forces[a_id];
                    Eigen::Vector3d torque = torques[a_id];
                    Eigen::Vector3d force_offset = force_offsets[a_id];
                    Eigen::Vector3d torque_offset = torque_offsets[a_id];
                    double t = s * cluster.step_size;
                    SStates[a_id] = updateState(hits[primary].p_a->current_state, t, force, force_offset, torque, torque_offset, hits[primary].p_a);
                }

                if (!hits[primary].p_b->is_static)
                {
                    Eigen::Vector3d force = forces[b_id];
                    Eigen::Vector3d torque = torques[b_id];
                    Eigen::Vector3d force_offset = force_offsets[b_id];
                    Eigen::Vector3d torque_offset = torque_offsets[b_id];
                    double t = s * cluster.step_size;
                    SStates[b_id] = updateState(hits[primary].p_b->current_state, t, force, force_offset, torque, torque_offset, hits[primary].p_b);
                }

                double outer = hits[primary].eps_a + hits[primary].eps_b;
                double inner = hits[primary].eps_inner_a + hits[primary].eps_inner_b;
                Eigen::Vector3d n = contacts[primary].global_normal.normalized();
                auto F = [=](double dist) { return std::max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))); };
                auto integral = [=](double dist) { return (inner - outer) / (dist - inner) - dist / (inner - outer); };

                Eigen::Vector3d force_last = -n * contacts[primary].force_last;
                Eigen::Vector3d pre_A = common::transform(contacts[primary].local_A, SStates[a_id].getTransformation());
                Eigen::Vector3d t_a_last = model::calcTorque(force_last, pre_A, cluster.particles[a_id].current_state.getTranslation());
                Eigen::Vector3d pre_B = common::transform(contacts[primary].local_B, SStates[b_id].getTransformation());
                Eigen::Vector3d t_b_last = model::calcTorque(Eigen::Vector3d(-force_last), pre_B, cluster.particles[b_id].current_state.getTranslation());

                Eigen::Vector3d p0_A = common::transform(contacts[primary].local_A, hits[primary].p_a->current_state.getTransformation());
                Eigen::Vector3d p0_B = common::transform(contacts[primary].local_B, hits[primary].p_b->current_state.getTransformation());
                double d0 = p0_B.dot(n) - p0_A.dot(n);

                Eigen::Vector3d ps_A = common::transform(contacts[primary].local_A, SStates[a_id].getTransformation());
                Eigen::Vector3d ps_B = common::transform(contacts[primary].local_B, SStates[b_id].getTransformation());
                double ds = ps_B.dot(n) - ps_A.dot(n);
                double ds_orig = ds;

                // printf("Start of iteration ds: %f\n", ds);

                Eigen::Vector3d v0_A = hits[primary].p_a->current_state.pointVelocity(p0_A, hits[primary].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v0_B = hits[primary].p_b->current_state.pointVelocity(p0_B, hits[primary].p_b->getInverseInertiaMatrix());
                double v0 = v0_A.dot(n) - v0_B.dot(n);

                Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[primary].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[primary].p_b->getInverseInertiaMatrix());
                double vs = vs_A.dot(n) - vs_B.dot(n);

                // printf("vs: %f\n", vs);

                if (v0 > 0.0) {
                    if (vs > 0.0) {
                        // Increase s
                        s += s_diff;
                    }
                    else if (vs < 0.0) {
                        // Decrease s
                        s -= s_diff;
                    }
                }

                s = std::min(std::max(s, 0.0), 1.0);

                s_diff *= s_diff_adjust;
                // printf("s_diff: %.10f\n", s_diff);

                contacts[primary].s = s;

                // printf("Picking S: %f\n", s);

                // double dt = cluster.step_size;
                // double F_ext = forces[b_id].dot(n) - forces[a_id].dot(n);
                // double Meff = contacts[primary].mass_eff_normal;
                // double F_last = contacts[primary].force_last;
                // // double F_last = 0.0;

                // auto Heaviside = [](double f) { return f >= 0.0 ? 1.0 : 0.0; };

                // auto eq = [=](double dist, double s) { return -s*(-v0 + vs)/s_last - v0 + dt*s*(-F_last + 0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2.0) - 1.0/(-inner + outer)) + 0.5*std::max(0.0, -1.0/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2.0)))/Meff; };
                // auto grad = [=](double dist, double s) { return -(-v0 + vs)/s_last - 1.0*dt*s*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(-F_last + 0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };

                // double no_force_intersection = s_last * (inner - d0) / (ds - d0);

                // if (no_force_intersection <= s && no_force_intersection > 0.0) {
                //     s = no_force_intersection * 0.9;
                //     double d = d0 + s * (ds - d0) / s_last;
                //     assert(grad(d, s) > 0.0);
                // }

                // printf("v0: %f, vs: %f, sl: %f, t: %f, Fl: %f, i: %f, o:%f, d0: %f, ds: %f, M: %f\n", v0, vs, s_last, dt, F_last, inner, outer, d0, ds, Meff);

                // for (int i = 0; i < s_calc_iterations; i++) {
                //     double d = d0 + s * (ds - d0) / s_last;
                //     double e = eq(d, s);
                //     double g = grad(d, s);
                //     s -= e / g;
                //     s = std::min(std::max(0.0, s), no_force_intersection > 0.0 && no_force_intersection < 1.0 ? no_force_intersection - 1e-3 : 1.0);
                // }

                // double d = d0 + s * (ds - d0) / s_last;
                // double e = eq(d, s);
                // double g = grad(d, s);
                // if (std::fabs(e) > 1e-6 && s < 1e-6 && (v0 > 0 || vs > 0)) {
                //     // There may still be a solution in the interval but we just started on the wrong side of the curve
                //     printf("Correcting s from: %f", s);
                //     s = no_force_intersection > 0.0 && no_force_intersection < 1.0 ? no_force_intersection * 0.99 : 1.0;

                //     for (int i = 0; i < s_calc_iterations; i++) {
                //         double d = d0 + s * (ds - d0) / s_last;
                //         double e = eq(d, s);
                //         double g = grad(d, s);
                //         s -= e / g;
                //         s = std::min(std::max(0.0, s), no_force_intersection > 0.0 && no_force_intersection < 1.0 ? no_force_intersection - 1e-3 : 1.0);
                //     }
                //     printf(" to %f\n", s);
                // }

                // printf("s: %f\n", s);
                // s = (1 - s_damp) * s_last + s_damp * s;
                // printf("s: %f (damped)\n", s);

                // Move all the particles to the new time step
                for (int p = 0; p < cluster.particles.size(); p++) {
                    if (!cluster.particles[p].is_static)
                    {
                        Eigen::Vector3d force = forces[p];
                        Eigen::Vector3d torque = torques[p];
                        Eigen::Vector3d force_offset = force_offsets[p];
                        Eigen::Vector3d torque_offset = torque_offsets[p];
                        double t = s * cluster.step_size;
                        // printf("Updating %i with t: %f, force: %f, %f, %f\n", p, t, force.x(), force.y(), force.z());
                        SStates[p] = updateState(cluster.particles[p].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[p]);
                    }
                }
            }

            for (int new_s_it = 0; new_s_it < new_S_iterations; new_s_it++) {
                for (int p = 0; p < cluster.particles.size(); p++) {
                    forces_last[p] = forces[p];
                    torques_last[p] = torques[p];
                    forces[p] = external_force(cluster.particles[p]);
                    torques[p] = { 0.0, 0.0, 0.0 };
                }

                // Update SStates for anything related to the contact we've just updated
                for (int c = 0; c < hits.size(); c++) {
                    uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                    uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                    double outer = hits[c].eps_a + hits[c].eps_b;
                    double inner = hits[c].eps_inner_a + hits[c].eps_inner_b;
                    Eigen::Vector3d n = contacts[c].global_normal.normalized();
                    auto F = [=](double dist) { return std::max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))); };
                    auto integral = [=](double dist) { return (inner - outer) / (dist - inner) - dist / (inner - outer); };

                    Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, hits[c].p_a->current_state.getTransformation());
                    Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, hits[c].p_b->current_state.getTransformation());
                    double d0 = p0_B.dot(n) - p0_A.dot(n);

                    Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, SStates[a_id].getTransformation());
                    Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, SStates[b_id].getTransformation());
                    double ds = ps_B.dot(n) - ps_A.dot(n);
                        
                    Eigen::Vector3d v0_A = hits[c].p_a->current_state.pointVelocity(p0_A, hits[c].p_a->getInverseInertiaMatrix());
                    Eigen::Vector3d v0_B = hits[c].p_b->current_state.pointVelocity(p0_B, hits[c].p_b->getInverseInertiaMatrix());
                    double v0 = v0_A.dot(n) - v0_B.dot(n);

                    Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
                    Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
                    double vs = vs_A.dot(n) - vs_B.dot(n);

                    double F0 = F(d0);
                    double F1 = F(ds);
                    contacts[c].force_last = averageForce(F0, F1, d0, ds, outer);

                    if (it == pre_S_iterations - 1 && new_s_it == new_S_iterations - 1) {
                        printf("Contact: %i, d0: %f, ds: %f, v0: %f, vs: %f\n", c, d0, ds, v0, vs);
                    }

                    if (v0 < 0 && vs < 0 && 1.0 - s < 1e-6) {
                        // printf("Damping from %f to %f (pre S)\n", contacts[c].force_last, forceDamp(contacts[c].force_last, v0, vs, separation_damp));
                        contacts[c].force_last = forceDamp(contacts[c].force_last, v0, vs, separation_damp);
                    }

                    Eigen::Vector3d force = -n * contacts[c].force_last;

                    // if (new_s_it == new_S_iterations - 1) {
                    //     printf("Contact: %i, force: %f, %f, %f, d0: %f, ds: %f\n", c, force.x(), force.y(), force.z(), d0, ds);
                    // }

                    Eigen::Vector3d T_A = model::calcTorque(force, ps_A, SStates[a_id].getTranslation());
                    Eigen::Vector3d T_B = model::calcTorque(Eigen::Vector3d(-force), ps_B, SStates[b_id].getTranslation());

                    forces[a_id] += force;
                    forces[b_id] -= force;
                    torques[a_id] += T_A;
                    torques[b_id] += T_B;
                }

                force_damp = 1.0 - double(new_s_it) / (double(new_S_iterations) - 1);

                for (int p = 0; p < cluster.particles.size(); p++) {
                    forces[p] = force_damp * forces[p] + (1.0 - force_damp) * forces_last[p];
                    torques[p] = force_damp * torques[p] + (1.0 - force_damp) * torques_last[p];
                }

                //force_damp *= force_damp_adjust;
                // printf("Force damp: %f\n", force_damp);

                for (int p = 0; p < cluster.particles.size(); p++) {
                    if (!cluster.particles[p].is_static)
                    {
                        Eigen::Vector3d force = forces[p];
                        Eigen::Vector3d torque = torques[p];
                        Eigen::Vector3d force_offset = force_offsets[p];
                        Eigen::Vector3d torque_offset = torque_offsets[p];
                        double t = s * cluster.step_size;
                        // printf("Updating end with t: %f\n", t);
                        SStates[p] = updateState(cluster.particles[p].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[p]);
                    }
                }
            }

            // select the next primary
            primary = 0;
            double next_primary_v = std::numeric_limits<double>::infinity();
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                if (!(a_id == primary_a_id || b_id == primary_a_id || a_id == primary_b_id || b_id == primary_b_id)) {
                    continue;
                }

                Eigen::Vector3d n = contacts[c].global_normal.normalized();

                Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, hits[c].p_a->current_state.getTransformation());
                Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, hits[c].p_b->current_state.getTransformation());

                Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, SStates[a_id].getTransformation());
                Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, SStates[b_id].getTransformation());

                Eigen::Vector3d v0_A = hits[c].p_a->current_state.pointVelocity(p0_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v0_B = hits[c].p_b->current_state.pointVelocity(p0_B, hits[c].p_b->getInverseInertiaMatrix());
                double v0 = v0_A.dot(n) - v0_B.dot(n);

                Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
                double vs = vs_A.dot(n) - vs_B.dot(n);

                if (v0 > 0.0) {
                    if (vs < next_primary_v) {
                        primary = c;
                        next_primary_v = vs;
                    }
                }
            }
            // printf("Selecting %i as the next primary with velocity %f\n", primary, next_primary_v);
        }
    }

    for (int c = 0; c < hits.size(); c++) {
        hits[c].force_mag = contacts[c].force_last;
        // printf("Normal: %f, %f, %f, Contact: %f, i\n", contacts[c].global_normal.x(), contacts[c].global_normal.y(), contacts[c].global_normal.z(), contacts[c].force_last, c);

        uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
        uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

        min_s[a_id] = std::min(min_s[a_id], contacts[c].s);
    }

    // Update all the states and make FState predictions ready for friction calculations
    std::vector<double> particleSStepSize;
    for (int p = 0; p < cluster.particles.size(); p++) {
        uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
        particleSStepSize.push_back(s * cluster.step_size);
        if (!cluster.particles[id].is_static) {
            Eigen::Vector3d force = forces[id];
            Eigen::Vector3d torque = torques[id];
            double t = s * cluster.step_size;
            Eigen::Vector3d force_offset = force_offsets[id];
            Eigen::Vector3d torque_offset = torque_offsets[id];
            printf("Particle: %i, SForce: %f, %f, %f, STorque: %f, %f, %f\n", p, force.x(), force.y(), force.z(), torque.x(), torque.y(), torque.z());
            State SState = updateState(cluster.particles[id].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            SStates[p] = SState;
            State FState = updateState(SState, (1 - s)*cluster.step_size, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            FStates[p] = FState;
            printf("0 vel: (%f, %f, %f), Z loc: %f\n", cluster.particles[id].current_state.getVelocity().x(), cluster.particles[id].current_state.getVelocity().y(), cluster.particles[id].current_state.getVelocity().z(), cluster.particles[id].current_state.getTranslation().z());
            printf("S vel: (%f, %f, %f), Z loc: %f, s: %f\n", SState.getVelocity().x(), SState.getVelocity().y(), SState.getVelocity().z(), SState.getTranslation().z(), s);
            printf("1 vel: (%f, %f, %f), Z loc: %f\n", FState.getVelocity().x(), FState.getVelocity().y(), FState.getVelocity().z(), FState.getTranslation().z());
            force_offsets[id] = {0, 0, 0};
            torque_offsets[id] = {0, 0, 0};
        }
        else {
            SStates[id] = cluster.particles[id].current_state;
            FStates[id] = cluster.particles[id].future_state;
        }
    }

    printf("Pre friction force:  %f, %f, %f, torque: %f, %f, %f\n", forces[0].x(), forces[0].y(), forces[0].z(), torques[0].x(), torques[0].y(), torques[0].z());
    _friction.solve(cluster.particles, hits, forces, torques, SStates, particleSStepSize);
    printf("Friction solve\n");
    printf("Post friction force: %f, %f, %f, torque: %f, %f, %f\n", forces[0].x(), forces[0].y(), forces[0].z(), torques[0].x(), torques[0].y(), torques[0].z());

    // Update states (including Fstate) ready for Fstate computation
    for (int p = 0; p < cluster.particles.size(); p++) {
        uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
        particleSStepSize.push_back(s * cluster.step_size);
        if (!cluster.particles[id].is_static) {
            Eigen::Vector3d force = forces[id];
            Eigen::Vector3d torque = torques[id];
            double t = s * cluster.step_size;
            Eigen::Vector3d force_offset = force_offsets[id];
            Eigen::Vector3d torque_offset = torque_offsets[id];
            State SState = updateState(cluster.particles[id].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            SStates[p] = SState;
            State FState = updateState(SState, (1.0 - s)*cluster.step_size, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            FStates[p] = FState;
            printf("0 vel: (%f, %f, %f), Z loc: %f\n", cluster.particles[id].current_state.getVelocity().x(), cluster.particles[id].current_state.getVelocity().y(), cluster.particles[id].current_state.getVelocity().z(), cluster.particles[id].current_state.getTranslation().z());
            printf("S vel: (%f, %f, %f), Z loc: %f, s: %f\n", SState.getVelocity().x(), SState.getVelocity().y(), SState.getVelocity().z(), SState.getTranslation().z(), s);
            printf("1 vel: (%f, %f, %f), Z loc: %f\n", FState.getVelocity().x(), FState.getVelocity().y(), FState.getVelocity().z(), FState.getTranslation().z());
            force_offsets[id] = {0, 0, 0};
            torque_offsets[id] = {0, 0, 0};
        }
        else {
            SStates[p] = cluster.particles[id].current_state;
            FStates[p] = cluster.particles[id].future_state;
        }
    }

    // Solve post S
    if (hits.size() > 0 && s < 1.0) {
        for (int it = 0; it < post_S_iterations; it++) {
            // for (int c = 0; c < hits.size(); c++) {
            //     uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            //     uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

            //     double outer = hits[c].eps_a + hits[c].eps_b;
            //     double inner = hits[c].eps_inner_a + hits[c].eps_inner_b;
            //     Eigen::Vector3d n = contacts[c].global_normal.normalized();
            //     auto F = [=](double dist) { return std::max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))); };

            //     Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, SStates[a_id].getTransformation());
            //     Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, SStates[b_id].getTransformation());
            //     // double ds = (ps_A - ps_B).norm();
            //     double ds = ps_A.dot(n) - ps_B.dot(n);

            //     Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            //     Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
            //     // double d1 = (p1_A - p1_B).norm();
            //     double d1 = p1_A.dot(n) - p1_B.dot(n);

            //     double dt = cluster.step_size;
            //     double Meff = contacts[c].mass_eff_normal;

            //     double F_last = contacts[c].force_last;

            //     Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
            //     Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
            //     double vs = -(vs_A.dot(n) - vs_B.dot(n));  // Negative compared to first step as optimisation uses +ve velocity to mean separating velocity

            //     double F_total = 0.0;  // TODO

            //     auto Heaviside = [](double f) { return f >= 0.0 ? 1.0 : 0.0; };
            //     // printf("ds: %f, dt: %f, vs: %f, s: %f, F_total: %f, inner: %f, outer: %f, Meff: %f\n", ds, dt, vs, s, F_total, inner, outer, Meff);
            //     auto eq = [=](double dc, double d1) { return d1 + dc - ds - dt*vs*(1 - s) - 0.5*std::pow(dt, 2)*std::pow(1 - s, 2)*(F_total + 0.5*std::max(0.0, (-inner + outer)/std::pow(-inner + std::max(0.0, std::min(ds, outer)), 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, (-inner + outer)/std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 2) - 1/(-inner + outer)))/Meff; };
            //     auto grad = [=](double dc) { return 1 + 0.5*std::pow(dt, 2)*std::pow(1 - s, 2)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 2) - 1/(-inner + outer))*Heaviside(-d1 - dc + outer)*Heaviside(std::min(outer, d1 + dc))/(Meff*std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 3)); };

            //     double d1_orig = d1;
            //     double dc = 0.0;
            //     for (int i = 0; i < post_S_distance_iterations; i++) {
            //         double e = eq(dc, d1);
            //         double g = grad(dc);
            //         // printf("dc: %f, d1: %f, e: %f, g: %f\n", dc, d1, e, g);
            //         dc -= e / g;
            //     }
            //     d1 = d1_orig + dc;

            //     printf("Contact: %i, ds: %f, d1_orig: %f, d1: %f, vs: %f\n", c, ds, d1_orig, d1, vs);

            //     // if (it == post_S_iterations - 1) {
            //     //     printf("Contact: %i, ds: %f, d1_orig: %f, d1: %f, vs: %f\n", c, ds, d1_orig, d1, vs);
            //     // }

            //     // printf("d1 estimate: %f\n", d1);

            //     double Fs = F(ds);
            //     double F1 = F(d1);

            //     // printf("ds: %f, d1: %f, Fs: %f, F1: %f, d1_orig: %f\n", ds, d1, Fs, F1, d1_orig);

            //     // printf("ds: %f, d1: %f, Fs: %f, F1: %f\n", ds, d1, Fs, F1);

            //     double last = contacts[c].force_last;
            //     contacts[c].force_last = averageForce(Fs, F1, ds, d1, outer);

            //     if (vs > 0 && 1.0 - s < 1e-6) {
            //         printf("Damping from %f to %f (post S)\n", contacts[c].force_last, contacts[c].force_last * separation_damp);
            //         contacts[c].force_last *= separation_damp;
            //     }

            //     // printf("Force: %f, last: %f\n", contacts[c].force_last, last);

            //     double force_diff = contacts[c].force_last - last;

            //     Eigen::Vector3d f = -n * force_diff;

            //     // printf("ps_A: %f, %f, %f\n", ps_A.x(), ps_A.y(), ps_A.z());

            //     Eigen::Vector3d t_a = model::calcTorque(f, ps_A, FStates[a_id].getTranslation());
            //     Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), ps_B, FStates[b_id].getTranslation());

            //     // printf("Torque: %f, %f, %f\n", t_a.x(), t_a.y(), t_a.z());

            //     forces[a_id] += f;
            //     torques[a_id] += t_a;
            //     forces[b_id] += -f;
            //     torques[b_id] += t_b;
            // }

            for (int p = 0; p < cluster.particles.size(); p++) {
                forces_last[p] = forces[p];
                torques_last[p] = torques[p];
                forces[p] = external_force(cluster.particles[p]);
                torques[p] = { 0.0, 0.0, 0.0 };
            }

            // Update SStates for anything related to the contact we've just updated
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                double outer = hits[c].eps_a + hits[c].eps_b;
                double inner = hits[c].eps_inner_a + hits[c].eps_inner_b;
                Eigen::Vector3d n = contacts[c].global_normal.normalized();
                auto F = [=](double dist) { return std::max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))); };
                auto integral = [=](double dist) { return (inner - outer) / (dist - inner) - dist / (inner - outer); };

                Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, SStates[a_id].getTransformation());
                Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, SStates[b_id].getTransformation());
                double ds = ps_B.dot(n) - ps_A.dot(n);
                    
                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
                double d1 = p1_B.dot(n) - p1_A.dot(n);

                Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
                double vs = vs_A.dot(n) - vs_B.dot(n);

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());
                double v1 = v1_A.dot(n) - v1_B.dot(n);

                double F0 = F(ds);
                double F1 = F(d1);
                contacts[c].force_last = averageForce(F0, F1, ds, d1, outer);

                if (it == post_S_iterations - 1) {
                    printf("Contact: %i, ds: %f, d1: %f, vs: %f, v1: %f\n", c, ds, d1, vs, v1);
                }

                if (v1 < 0 && 1.0 - s > 1e-6) {
                    // if (it == post_S_iterations - 1) {
                    //     printf("Damping from %f to %f (post S)\n", contacts[c].force_last, forceDamp(contacts[c].force_last, vs, v1, separation_damp));
                    // }
                    contacts[c].force_last = forceDamp(contacts[c].force_last, vs, v1, separation_damp);
                    // contacts[c].force_last *= separation_damp;
                }

                Eigen::Vector3d force = -n * contacts[c].force_last;

                // if (new_s_it == new_S_iterations - 1) {
                //     printf("Contact: %i, force: %f, %f, %f, d0: %f, ds: %f\n", c, force.x(), force.y(), force.z(), d0, ds);
                // }

                Eigen::Vector3d T_A = model::calcTorque(force, p1_A, FStates[a_id].getTranslation());
                Eigen::Vector3d T_B = model::calcTorque(Eigen::Vector3d(-force), p1_B, FStates[b_id].getTranslation());

                forces[a_id] += force;
                forces[b_id] -= force;
                torques[a_id] += T_A;
                torques[b_id] += T_B;
            }

            force_damp = 1.0 - double(it) / (double(post_S_iterations) - 1);

            for (int p = 0; p < cluster.particles.size(); p++) {
                forces[p] = force_damp * forces[p] + (1.0 - force_damp) * forces_last[p];
                torques[p] = force_damp * torques[p] + (1.0 - force_damp) * torques_last[p];
            }

            for (int p = 0; p < cluster.particles.size(); p++) {
                uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
                if (!cluster.particles[id].is_static) {
                    Eigen::Vector3d force = forces[id];
                    Eigen::Vector3d torque = torques[id];
                    double t = (1.0 - s) * cluster.step_size;
                    Eigen::Vector3d force_offset = { 0.0, 0.0, 0.0 };
                    Eigen::Vector3d torque_offset = { 0.0, 0.0, 0.0 };

                    FStates[id] = updateState(SStates[id], t, force, force_offset, torque, torque_offset, &cluster.particles[p]);
                }
            }
        }

        printf("Post S contacts\n");
        for (int c = 0; c < hits.size(); c++) {
            printf("Contact: %i, Normal: %f, %f, %f, Force: %f\n", c, contacts[c].global_normal.x(), contacts[c].global_normal.y(), contacts[c].global_normal.z(), contacts[c].force_last);
        }

        // for (int p = 0; p < cluster.particles.size(); p++) {
        //     forces[p] = external_force(cluster.particles[p]);
        //     torques[p] = { 0.0, 0.0, 0.0 };
        // }

        // for (int c = 0; c < hits.size(); c++) {
        //     uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
        //     uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

        //     Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, SStates[a_id].getTransformation());
        //     Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, SStates[b_id].getTransformation());

        //     Eigen::Vector3d f = -contacts[c].global_normal * contacts[c].force_last * separation_damp;

        //     Eigen::Vector3d t_a = model::calcTorque(f, ps_A, FStates[a_id].getTranslation());
        //     Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), ps_B, FStates[b_id].getTranslation());

        //     forces[a_id] += f;
        //     torques[a_id] += t_a;
        //     forces[b_id] += -f;
        //     torques[b_id] += t_b;
        // }

        // for (int p = 0; p < cluster.particles.size(); p++) {
        //     uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
        //     if (!cluster.particles[id].is_static) {
        //         Eigen::Vector3d force = forces[id];
        //         Eigen::Vector3d torque = torques[id];
        //         double t = (1.0 - min_s[id]) * cluster.step_size;
        //         Eigen::Vector3d force_offset = { 0.0, 0.0, 0.0 };
        //         Eigen::Vector3d torque_offset = { 0.0, 0.0, 0.0 };

        //         FStates[id] = updateState(SStates[id], t, force, force_offset, torque, torque_offset, &cluster.particles[p]);
                
        //         // FStates[id].setVelocity(FStates[id].getVelocity() * separation_damp);
        //         // FStates[id].setAngular(FStates[id].getAngularMomentum() * separation_damp);
        //     }
        // }


        _friction.solve(cluster.particles, hits, forces, torques, FStates, particleSStepSize);

        // Update states (including Fstate) ready for Fstate computation
        for (int p = 0; p < cluster.particles.size(); p++) {
            uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
            particleSStepSize.push_back(s * cluster.step_size);
            if (!cluster.particles[id].is_static) {
                Eigen::Vector3d force = forces[id];
                Eigen::Vector3d torque = torques[id];
                Eigen::Vector3d force_offset = force_offsets[id];
                Eigen::Vector3d torque_offset = torque_offsets[id];
                State FState = updateState(SStates[id], (1.0 - s)*cluster.step_size, force, force_offset, torque, torque_offset, &cluster.particles[id]);
                FStates[p] = FState;
            }
            else {
                FStates[p] = cluster.particles[id].future_state;
            }
        }
    }

    double drag = std::pow(0.5, cluster.step_size);
    drag = 1;

    for (Particle *p : cluster.particles)
    {
        if (!p->is_static) {
            uint32_t id = cluster.particles.getLocalID(p->id);
            FStates[id].setAngular(FStates[id].getAngularMomentum() * drag);
            FStates[id].setVelocity(FStates[id].getVelocity() * drag);
            p->future_state = FStates[id];
            printf("End of step z vel: %f, z loc: %f\n", p->future_state.getVelocity().z(), p->future_state.getTranslation().z());
        }
    }
    
    // common::Viewer view;
    // view.addParticleInterval(cluster.particles[0]);
    // view.addParticleInterval(cluster.particles[1]);
    // for (int c = 0; c < hits.size(); c++) {
    //     view.addEdge(common::Edge(hits[c].A, hits[c].B));
    // }
    // view.show();
}

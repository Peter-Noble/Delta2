#include "sequential_impulses_s_time.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"

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

    double force_last;
    double force_move_offset;
};

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

void SequentialImpulsesSTime::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) {
    std::vector<LocalContactS> contacts;

    std::sort(hits.begin(), hits.end());

    for (int c = 0; c < hits.size(); c++) {
        LocalContactS lc;

        Eigen::Matrix4d A_T_i = hits[c].p_a->future_state.getTransformation().inverse();
        lc.local_A = common::transform(hits[c].A, A_T_i);
        Eigen::Matrix4d B_T_i = hits[c].p_b->future_state.getTransformation().inverse();
        lc.local_B = common::transform(hits[c].B, B_T_i);

        lc.global_centre = (hits[c].A + hits[c].B) / 2.0;
        lc.global_normal = (hits[c].B - hits[c].A) / 2.0;
        lc.global_tangent = {0, 0, 0}; // TODO

        printf("Normal: %f, %f, %f\n", lc.global_normal.x(), lc.global_normal.y(), lc.global_normal.z());

        if (lc.global_normal.norm() < 1e-6) {
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

    std::vector<Eigen::Vector3d> impulses;
    std::vector<Eigen::Vector3d> rotational_impulses;
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<Eigen::Vector3d> force_offsets;
    std::vector<Eigen::Vector3d> torque_offsets;
    std::vector<double> min_s;
    min_s.resize(cluster.particles.size());

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        forces.push_back(external_force(cluster.particles[p_i]));
        torques.push_back({0, 0, 0});
        force_offsets.push_back({0, 0, 0});
        torque_offsets.push_back({0, 0, 0});
    }

    double last_ds;

    for (int it = 0; it < 10; it++) {
        for (int p = 0; p < cluster.particles.size(); p++) {
            min_s[p] = 1.0;
        }
        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

            // Integrate without the contribution of this contact.  Then solve for s, work out the contribution of this contact and add it.

            double outer = hits[c].eps_a + hits[c].eps_b;
            double inner = hits[c].eps_inner_a + hits[c].eps_inner_b;
            Eigen::Vector3d n = contacts[c].global_normal.normalized();
            auto F = [=](double dist) { return std::max(0.0, (outer - inner) / ((dist - inner)*(dist - inner)) - (outer - inner) / ((outer - inner)*(outer - inner))); };
            auto integral = [=](double dist) { return (inner - outer) / (dist - inner) - dist / (inner - outer); };

            Eigen::Vector3d force_last = -n * contacts[c].force_last;
            Eigen::Vector3d pre_A = common::transform(contacts[c].local_A, contacts[c].S_A.getTransformation());
            Eigen::Vector3d t_a_last = model::calcTorque(force_last, pre_A, cluster.particles[a_id].current_state.getTranslation());
            Eigen::Vector3d pre_B = common::transform(contacts[c].local_B, contacts[c].S_B.getTransformation());
            Eigen::Vector3d t_b_last = model::calcTorque(Eigen::Vector3d(-force_last), pre_B, cluster.particles[b_id].current_state.getTranslation());

            if (!hits[c].p_a->is_static)
            {
                // Eigen::Vector3d force = forces[a_id] - force_last;
                // Eigen::Vector3d torque = torques[a_id] - t_a_last;
                Eigen::Vector3d force = forces[a_id]; // - force_last;
                Eigen::Vector3d torque = torques[a_id]; // - t_a_last;
                Eigen::Vector3d force_offset = force_offsets[a_id];
                Eigen::Vector3d torque_offset = torque_offsets[a_id];
                double t = contacts[c].s * cluster.step_size;
                contacts[c].S_A = updateState(hits[c].p_a->current_state, t, force, force_offset, torque, torque_offset, hits[c].p_a);
                // contacts[c].S_A = updateState(hits[c].p_a->current_state, cluster.step_size, force, force_offset, torque, torque_offset, hits[c].p_a);
            }

            if (!hits[c].p_b->is_static)
            {
                // Eigen::Vector3d force = forces[b_id] + force_last;
                // Eigen::Vector3d torque = torques[b_id] - t_b_last;
                Eigen::Vector3d force = forces[b_id]; // + force_last;
                Eigen::Vector3d torque = torques[b_id]; // - t_b_last;
                Eigen::Vector3d force_offset = force_offsets[b_id];
                Eigen::Vector3d torque_offset = torque_offsets[b_id];
                double t = contacts[c].s * cluster.step_size;
                contacts[c].S_B = updateState(hits[c].p_b->current_state, t, force, force_offset, torque, torque_offset, hits[c].p_b);
                // contacts[c].S_B = updateState(hits[c].p_b->current_state, cluster.step_size, force, force_offset, torque, torque_offset, hits[c].p_b);
            }

            Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, hits[c].p_a->current_state.getTransformation());
            Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, hits[c].p_b->current_state.getTransformation());
            double d0 = p0_B.dot(n) - p0_A.dot(n);

            Eigen::Vector3d ps_A = common::transform(contacts[c].local_A, contacts[c].S_A.getTransformation());
            Eigen::Vector3d ps_B = common::transform(contacts[c].local_B, contacts[c].S_B.getTransformation());
            double ds = ps_B.dot(n) - ps_A.dot(n);
            double ds_orig = ds;

            // printf("Start of iteration ds: %f\n", ds);

            Eigen::Vector3d v0_A = hits[c].p_a->current_state.pointVelocity(p0_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d v0_B = hits[c].p_b->current_state.pointVelocity(p0_B, hits[c].p_b->getInverseInertiaMatrix());
            double v0 = v0_A.dot(n) - v0_B.dot(n);

            Eigen::Vector3d vs_A = contacts[c].S_A.pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d vs_B = contacts[c].S_B.pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
            double vs = vs_A.dot(n) - vs_B.dot(n);

            double dt = cluster.step_size;
            double F_ext = forces[b_id].dot(n) - forces[a_id].dot(n);
            double Meff = contacts[c].mass_eff_normal;
            double F_last = contacts[c].force_last;
            // double F_last = 0.0;

            double s_last = contacts[c].s;
            // s_last = 1.0;
            double s = s_last;

            auto Heaviside = [](double f) { return f >= 0.0 ? 1.0 : 0.0; };
            // These ones solve for final velocity of vs but the velocity changes over the time interval so the velocity we solve for needs to depend on s
            // auto eq = [=](double dist, double s) { return -vs + s*dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };
            // auto grad = [=](double dist, double s) { return -1.0*s*dt*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };

            // Works for the first iteration but then the velocity is zero so contact moves earlier and earlier
            // auto eq = [=](double dist, double s) { return -s*(-v0 + vs)/s_last - v0 + s*dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };
            // auto grad = [=](double dist, double s) { return -(-v0 + vs)/s_last - 1.0*s*dt*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };
 
            // auto eq = [=](double dist, double s) { return (-F_last*s_last*dt/Meff) - s*(-v0 + vs)/s_last - v0 + dt*s*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };
            // auto grad = [=](double dist, double s) { return (-F_last*dt/Meff) - (-v0 + vs)/s_last - 1.0*dt*s*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };

            // auto eq = [=](double dist, double s) { return -(s*(-F_last*dt*s/Meff - v0 + vs) + v0)/s_last + dt*s*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };
            // auto grad = [=](double dist, double s) { return -(-2*F_last*dt*s/Meff - v0 + vs)/s_last - 1.0*dt*s*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };

            auto eq = [=](double dist, double s) { return -s*(-v0 + vs)/s_last - v0 + dt*s*(-F_last + 0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2.0) - 1.0/(-inner + outer)) + 0.5*std::max(0.0, -1.0/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2.0)))/Meff; };
            auto grad = [=](double dist, double s) { return -(-v0 + vs)/s_last - 1.0*dt*s*(-d0 + ds)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer))/(Meff*s_last*std::pow(d0 - inner + s*(-d0 + ds)/s_last, 3)) + dt*(-F_last + 0.5*std::max(0.0, (-inner + outer)/std::pow(d0 - inner + s*(-d0 + ds)/s_last, 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, -1/(-inner + outer) + (-inner + outer)/std::pow(d0 - inner, 2)))/Meff; };

            double no_force_intersection = s_last * (inner - d0) / (ds - d0);

            if (no_force_intersection <= s) {
                s = no_force_intersection * 0.9;
                double d = d0 + s * (ds - d0) / s_last;
                assert(grad(d, s) > 0.0);
            }

            for (int i = 0; i < 10; i++) {
                double velocity_change = ((F(d0+s*(ds-d0)/s_last) + F(d0))/2.0 - F_last) * (s * dt) / Meff;
                double eq_rhs = (v0+s*(vs-v0)/s_last);

                double d = d0 + s * (ds - d0) / s_last;
                // printf("d: %f\n", d);
                double e = eq(d, s);
                // double e = velocity_change_sympy - eq_rhs_sympy;
                double g = grad(d, s);
                // printf("s: %f, e: %f, g %f\n", s, e, g);
                s -= e / g;
                s = std::min(std::max(0.0, s), 1.0);
                double force_z = (F(d0) + F(d)) / 2.0;
                // printf("Force z: %f\n", force_z);
                // printf("Eq: %f, %f\n", velocity_change, eq_rhs);
            }

            double last = contacts[c].force_last;

            if (s > 0.0) {
                contacts[c].s = s;
                min_s[a_id] = std::min(min_s[a_id], s);
                min_s[b_id] = std::min(min_s[b_id], s);

                ds = d0+s*(ds-d0)/s_last;
                // printf("New ds: %f\n", ds);
                double F0 = F(d0);
                double F1 = F(ds);
                contacts[c].force_last = std::max(0.0, (F(d0) + F(ds)) / 2.0);
            }
            else {
                contacts[c].s = 1.0;
                ds = d0+1.0*(ds-d0)/s_last;
                contacts[c].force_last = std::max(0.0, (F(d0) + F(ds)) / 2.0);
            }

            double force_diff = (contacts[c].force_last - last) * 0.8;
            contacts[c].force_last -= (contacts[c].force_last - last) * 0.2;

            Eigen::Vector3d f = -n * force_diff;

            Eigen::Vector3d t_a = model::calcTorque(f, ps_A, contacts[c].S_A.getTranslation());
            Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), ps_B, contacts[c].S_B.getTranslation());

            forces[a_id] += f;
            // printf("Force A: %f\n", forces[a_id].z());
            torques[a_id] += t_a;
            forces[b_id] += -f;
            torques[b_id] += t_b;

            Vector3d f_offset = -n * contacts[c].force_move_offset;
            Vector3d t_offset_a = model::calcTorque(f_offset, ps_A, contacts[c].S_A.getTranslation());  // TODO last two args here could probably be updated
            Vector3d t_offset_b = model::calcTorque(Vector3d(-f_offset), ps_B, contacts[c].S_B.getTranslation());  // TODO last two args here could probably be updated
            
            force_offsets[a_id] -= f_offset;
            force_offsets[b_id] += f_offset;
            torque_offsets[a_id] -= t_offset_a;
            torque_offsets[b_id] -= t_offset_b;

            if (s > 0) {
                contacts[c].force_move_offset = 0;
            }
            else {
                contacts[c].force_move_offset = 0;
            }

            f_offset = -n * contacts[c].force_move_offset;
            t_offset_a = model::calcTorque(f_offset, ps_A, contacts[c].S_A.getTranslation());  // TODO last two args here could probably be updated
            t_offset_b = model::calcTorque(Vector3d(-f_offset), ps_B, contacts[c].S_B.getTranslation());  // TODO last two args here could probably be updated
            
            // TODO Has this correction factor been correctly removed?  The two step approach will take care of this instead?  Or is this to take into account the changing velocity?
            force_offsets[a_id] += f_offset;
            force_offsets[b_id] -= f_offset;
            torque_offsets[a_id] += t_offset_a;
            torque_offsets[b_id] += t_offset_b;

            last_ds = ds;

            {
                Eigen::Vector3d force = forces[a_id];
                Eigen::Vector3d torque = torques[a_id];
                double t = min_s[a_id] * cluster.step_size;
                Eigen::Vector3d force_offset = force_offsets[a_id];
                Eigen::Vector3d torque_offset = torque_offsets[a_id];
                State SState = updateState(cluster.particles[a_id].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[a_id]);
                if (s < 0.9 && s > 0.1 && SState.getVelocity().norm() > 1e-2 && it == 9 ) {
                    printf("Strange\n");
                }
            }
        }
    }

    std::vector<State> SStates;
    std::vector<State> FStates;
    for (int p = 0; p < cluster.particles.size(); p++) {
        uint32_t id = cluster.particles.getLocalID(cluster.particles[p].id);
        if (!cluster.particles[id].is_static) {
            Eigen::Vector3d force = forces[id];
            Eigen::Vector3d torque = torques[id];
            double t = min_s[id] * cluster.step_size;
            Eigen::Vector3d force_offset = force_offsets[id];
            Eigen::Vector3d torque_offset = torque_offsets[id];
            State SState = updateState(cluster.particles[id].current_state, t, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            SStates.push_back(SState);
            State FState = updateState(SState, (1 - min_s[id])*cluster.step_size, force, force_offset, torque, torque_offset, &cluster.particles[id]);
            FStates.push_back(FState);
            printf("0 vel: (%f, %f, %f), Z loc: %f\n", cluster.particles[id].current_state.getVelocity().x(), cluster.particles[id].current_state.getVelocity().y(), cluster.particles[id].current_state.getVelocity().z(), cluster.particles[id].current_state.getTranslation().z());
            printf("S vel: (%f, %f, %f), Z loc: %f, s: %f\n", SState.getVelocity().x(), SState.getVelocity().y(), SState.getVelocity().z(), SState.getTranslation().z(), min_s[id]);
            printf("1 vel: (%f, %f, %f), Z loc: %f\n", FState.getVelocity().x(), FState.getVelocity().y(), FState.getVelocity().z(), FState.getTranslation().z());
            force_offsets[id] = {0, 0, 0};
            torque_offsets[id] = {0, 0, 0};
        }
        else {
            SStates.push_back(cluster.particles[id].current_state);
            FStates.push_back(cluster.particles[id].future_state);
        }
    }

    _friction.solve(ph, hits, forces, torques, counts, external_force, step_size);
    // TODO solve friction

    // Solve post S
    for (int it = 0; it < 10; it++) {
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
            double ds = (ps_A - ps_B).norm();

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
            double d1 = (p1_A - p1_B).norm();

            double dt = cluster.step_size;
            double Meff = contacts[c].mass_eff_normal;

            double F_last = contacts[c].force_last;
            double s = contacts[c].s;

            Eigen::Vector3d vs_A = SStates[a_id].pointVelocity(ps_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d vs_B = SStates[b_id].pointVelocity(ps_B, hits[c].p_b->getInverseInertiaMatrix());
            double vs = -(vs_A.dot(n) - vs_B.dot(n));  // Negative compared to first step as optimisation uses +ve velocity to mean separating velocity

            double F_total = 0.0;  // TODO

            auto Heaviside = [](double f) { return f >= 0.0 ? 1.0 : 0.0; };
            auto eq = [=](double dc, double d1) { return d1 + dc - ds - dt*vs*(1 - s) - 0.5*std::pow(dt, 2)*std::pow(1 - s, 2)*(F_total + 0.5*std::max(0.0, (-inner + outer)/std::pow(-inner + std::max(0.0, std::min(ds, outer)), 2) - 1/(-inner + outer)) + 0.5*std::max(0.0, (-inner + outer)/std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 2) - 1/(-inner + outer)))/Meff; };
            auto grad = [=](double dc) { return 1 + 0.5*std::pow(dt, 2)*std::pow(1 - s, 2)*(-inner + outer)*Heaviside((-inner + outer)/std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 2) - 1/(-inner + outer))*Heaviside(-d1 - dc + outer)*Heaviside(std::min(outer, d1 + dc))/(Meff*std::pow(-inner + std::max(0.0, std::min(outer, d1 + dc)), 3)); };

            double d1_orig = d1;
            double dc = 0.0;
            for (int i = 0; i < 10; i++) {
                double e = eq(dc, d1);
                double g = grad(dc);
                dc -= e / g;
            }
            d1 = d1_orig + dc;

            // printf("d1 estimate: %f\n", d1);

            double Fs = F(ds);
            double F1 = F(d1);

            double last = contacts[c].force_last;
            contacts[c].force_last = std::max(0.0, ((Fs + F1) / 2.0));

            double force_diff = contacts[c].force_last - last;

            Eigen::Vector3d f = -n * force_diff;

            Eigen::Vector3d t_a = model::calcTorque(f, ps_A, contacts[c].S_A.getTranslation());
            Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), ps_B, contacts[c].S_B.getTranslation());

            forces[a_id] += f;
            torques[a_id] += t_a;
            forces[b_id] += -f;
            torques[b_id] += t_b;

            double t = (1-contacts[c].s) * cluster.step_size;

            if (!hits[c].p_a->is_static)
            {
                Eigen::Vector3d force = forces[a_id];
                Eigen::Vector3d torque = torques[a_id];
                Eigen::Vector3d force_offset = force_offsets[a_id];
                Eigen::Vector3d torque_offset = torque_offsets[a_id];
                FStates[a_id] = updateState(SStates[a_id], t, force, force_offset, torque, torque_offset, hits[c].p_a);
            }

            if (!hits[c].p_b->is_static)
            {
                Eigen::Vector3d force = forces[b_id];
                Eigen::Vector3d torque = torques[b_id];
                Eigen::Vector3d force_offset = force_offsets[b_id];
                Eigen::Vector3d torque_offset = torque_offsets[b_id];
                FStates[b_id] = updateState(SStates[b_id], t, force, force_offset, torque, torque_offset, hits[c].p_b);
            }
        }
    }


    for (Particle *p : cluster.particles)
    {
        if (!p->is_static) {
            uint32_t id = cluster.particles.getLocalID(p->id);
            p->future_state = FStates[id];
            printf("End of step z vel: %f, z loc: %f\n", p->future_state.getVelocity().z(), p->future_state.getTranslation().z());
        }

        // if (!p->is_static)
        // {
        //     uint32_t id = cluster.particles.getLocalID(p->id);
        //     Eigen::Vector3d force = forces[id];
        //     Eigen::Vector3d torque = torques[id];
        //     double diff = model::integrateEuler(*p, force + external_force(*p), torque, cluster.step_size, false);
        //     if (diff < 2.0 * p->geo_eps) {
        //         if (p->current_state.getTime() - p->sleep_candidate_time > 0.25) {
        //             p->setSleeping(true);
        //             printf("Particle: %i put to sleep\n", p->id);
        //         }
        //     }
        //     else {
        //         if (p->getSleeping()) {
        //             p->setSleeping(false);
        //             printf("Particle: %i woken up\n", p->id);
        //         }
        //         p->sleep_candidate_time = p->current_state.getTime();
        //     }
        // }
    }
}

#include "sequential_impulses.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"
#include "../../model/friction.h"

#include "../../common/viewer.h"

#include <algorithm>

using namespace Delta2;
using namespace strategy;
using namespace Eigen;

struct LocalContact {
    Eigen::Vector3d local_A;
    Eigen::Vector3d local_B;
    Eigen::Vector3d global_centre;
    Eigen::Vector3d global_normal;
    Eigen::Vector3d global_tangent;
    double mass_eff_normal;
    double mass_eff_tangent;

    double start_dist;

    double impulse;
    double impulse_offset;

    Eigen::Vector3d last_friction_impulse;
    Eigen::Vector3d last_friction_rotational_impulse_A;
    Eigen::Vector3d last_friction_rotational_impulse_B;
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

State updateState(State& future, double t, const Vector3d& impulse, const Vector3d& impulse_offset, const Vector3d& rotational_impulse, const Vector3d& rotational_impulse_offset, Particle* p) {
    State updated = future;

    updated.setVelocity(future.getVelocity() + impulse / p->getMass());
    updated.setTranslation(future.getTranslation() + impulse_offset / p->getMass() * t);

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = future.getAngularMomentum() + rotational_impulse;

    Eigen::Quaterniond rotation = future.getRotation();

    updated.setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * p->getInverseInertiaMatrix() * R.transpose();
    Eigen::Vector3d omega = Iinv * (future.getAngularMomentum() + t * rotational_impulse_offset);

    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;
    updated.setRotation(rotation);

    return updated;
}

SequentialImpulses::SequentialImpulses(FrictionStrategy& friction,
                                     common::Options& opt) :
                                     ContactForceStrategy(friction, opt) {

}

bool SequentialImpulses::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) {
    std::vector<LocalContact> contacts;

    std::sort(hits.begin(), hits.end());

    for (int c = 0; c < hits.size(); c++) {
        LocalContact lc;

        lc.start_dist = (hits[c].A - hits[c].B).norm();

        Eigen::Matrix4d A_T_i = hits[c].p_a->future_state.getTransformation().inverse();
        lc.local_A = common::transform(hits[c].A, A_T_i);
        Eigen::Matrix4d B_T_i = hits[c].p_b->future_state.getTransformation().inverse();
        lc.local_B = common::transform(hits[c].B, B_T_i);

        lc.global_centre = (hits[c].A + hits[c].B) / 2.0;
        lc.global_normal = (hits[c].B - hits[c].A) / 2.0;
        lc.global_tangent = {0, 0, 0}; // todo

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

        lc.impulse = 0.0;
        lc.impulse_offset = 0.0;

        lc.last_friction_impulse = { 0.0, 0.0, 0.0 };
        lc.last_friction_rotational_impulse_A = { 0.0, 0.0, 0.0 };

        contacts.push_back(lc);
    }

    std::vector<Eigen::Vector3d> impulses;
    std::vector<Eigen::Vector3d> rotational_impulses;
    std::vector<Eigen::Vector3d> impulse_offsets;
    std::vector<Eigen::Vector3d> rotational_impulse_offsets;

    std::vector<State> FStates;

    const double t = cluster.step_size;

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        impulses.push_back(external_force(cluster.particles[p_i]) * t);
        // impulses.push_back({0, 0, 0});
        rotational_impulses.push_back({0, 0, 0});
        impulse_offsets.push_back({0, 0, 0});
        rotational_impulse_offsets.push_back({0, 0, 0});
        FStates.push_back(cluster.particles[p_i].future_state);
        
        if (!cluster.particles[p_i].is_static)
        {
            Eigen::Vector3d impulse = impulses[p_i];
            Eigen::Vector3d rotational_impulse = rotational_impulses[p_i];
            Eigen::Vector3d impulse_offset = impulse_offsets[p_i];
            Eigen::Vector3d rotational_impulse_offset = rotational_impulse_offsets[p_i];
            FStates[p_i] = updateState(cluster.particles[p_i].future_state, t, impulse, impulse_offset, rotational_impulse, rotational_impulse_offset, &cluster.particles[p_i]);
        }
    }

    if (hits.size() > 0) {
        bool converged = false;
        int it = 0;
        while (!converged && it < 1000) {
            std::vector<Eigen::Vector3d> last_impulse;
            for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                last_impulse.push_back(impulses[p_i]);
            }
            // printf("it: %i\n", it);
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                Eigen::Vector3d n = contacts[c].global_normal.normalized();

                const double outer_search = hits[c].eps_a + hits[c].eps_b;
                const double outer = outer_search * 0.1;

                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
                const double d1 = p1_B.dot(n) - p1_A.dot(n);

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());
                const double v1 = v1_A.dot(n) - v1_B.dot(n);

                const double m = contacts[c].mass_eff_normal;

                // ======== CANCEL RELATIVE VELOCITY ALONG NORMAL ========

                const double delta_damp = 0.1;
                // v1 = im / m to cancel velocity
                const double delta_impulse_mag = d1 < outer ? v1 * m * delta_damp : -contacts[c].impulse;

                double last_impulse_mag = contacts[c].impulse;
                contacts[c].impulse = std::max(contacts[c].impulse + delta_impulse_mag, 0.0);
                const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;

                Eigen::Vector3d update_impulse = -n * update_impulse_mag;

                impulses[a_id] += update_impulse;
                impulses[b_id] -= update_impulse;

                Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, FStates[a_id].getTranslation());
                Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, FStates[b_id].getTranslation());

                rotational_impulses[a_id] += update_rotational_impulse_A;
                rotational_impulses[b_id] += update_rotational_impulse_B;

                // ========  OFFSET ALONG NORMAL ========

                const double lower_slop = std::max(outer * 0.1, outer - contacts[c].start_dist);
                const double slop = common::lerp(lower_slop, outer * 0.1, 0.3);
                const double penetration_depth = outer - d1;

                // d1 = i/m * t offset distance with impulse.  m * d1 / t = i
                const double delta_impulse_offset_mag = d1 < outer ? delta_damp * m * std::max(0.0, penetration_depth - slop) / t : -contacts[c].impulse_offset;

                double last_impulse_offset_mag = contacts[c].impulse_offset;
                contacts[c].impulse_offset = std::max(contacts[c].impulse_offset + delta_impulse_offset_mag, 0.0);
                const double update_impulse_offset_mag = contacts[c].impulse_offset - last_impulse_offset_mag;

                Eigen::Vector3d update_impulse_offset = -n * update_impulse_offset_mag;

                impulse_offsets[a_id] += update_impulse_offset;
                impulse_offsets[b_id] -= update_impulse_offset;

                Eigen::Vector3d update_rotational_impulse_offset_A = model::calcTorque(update_impulse_offset, p1_A, FStates[a_id].getTranslation());
                Eigen::Vector3d update_rotational_impulse_offset_B = model::calcTorque(Eigen::Vector3d(-update_impulse_offset), p1_B, FStates[b_id].getTranslation());

                rotational_impulse_offsets[a_id] += update_rotational_impulse_offset_A;
                rotational_impulse_offsets[b_id] += update_rotational_impulse_offset_B;

                // printf("contact: %i, v1: %f, impulse: %f, d1: %f, impulse offset: %f\n", c, v1, contacts[c].impulse, d1, contacts[c].impulse_offset);

                // ======== FRICTION ========

                Eigen::Vector3d v1_AB = v1_A - v1_B;

                Eigen::Vector3d v_tangent = v1_AB - v1_AB.dot(n) * n;
                double mu = sqrt(cluster.particles[a_id].friction_coeff * cluster.particles[a_id].friction_coeff + cluster.particles[b_id].friction_coeff * cluster.particles[b_id].friction_coeff);
                double max_force = contacts[c].impulse * mu;

                Eigen::Vector3d r_A = hits[c].A - hits[c].p_a->future_state.getTranslation();
                Eigen::Vector3d r_B = hits[c].B - hits[c].p_b->future_state.getTranslation();

                Eigen::Vector3d rn_A = r_A.cross(contacts[c].global_normal);
                Eigen::Vector3d rn_B = r_B.cross(contacts[c].global_normal);

                Eigen::Vector3d tangent = v_tangent.normalized();

                Eigen::Vector3d rt_A = r_A.cross(tangent);
                Eigen::Vector3d rt_B = r_B.cross(tangent);

                double im_A = 1.0 / hits[c].p_a->getMass();
                double im_B = 1.0 / hits[c].p_b->getMass();
                Eigen::Matrix3d ii_A = hits[c].p_a->getInverseInertiaMatrix();
                Eigen::Matrix3d ii_B = hits[c].p_b->getInverseInertiaMatrix();

                double k_tangent = 0.0;
                if (!hits[c].p_a->is_static) {
                    k_tangent += im_A + rt_A.transpose() * ii_A * rt_A;

                }
                if (!hits[c].p_b->is_static) {
                    k_tangent += im_B + rt_B.transpose() * ii_B * rt_B;
                }

                double m_inv_effective_tangent = k_tangent < 1e-6 ? 0.0 : 1.0 / k_tangent;

                Eigen::Vector3d delta_friction_impulse = -v_tangent * m_inv_effective_tangent;
                Eigen::Vector3d friction_impulse_total = contacts[c].last_friction_impulse + delta_friction_impulse;
                friction_impulse_total = friction_impulse_total.normalized() * std::clamp(friction_impulse_total.norm(), 0.0, max_force);

                Eigen::Vector3d update_friction_impulse = friction_impulse_total - contacts[c].last_friction_impulse;

                contacts[c].last_friction_impulse = friction_impulse_total;
                impulses[a_id] += update_friction_impulse;
                impulses[b_id] -= update_friction_impulse;

                Eigen::Vector3d rotational_impulse_friction_A = model::calcTorque(contacts[c].last_friction_impulse, p1_A, FStates[a_id].getTranslation());
                Eigen::Vector3d rotational_impulse_friction_B = model::calcTorque(Eigen::Vector3d(-contacts[c].last_friction_impulse), p1_B, FStates[b_id].getTranslation());

                rotational_impulses[a_id] -= contacts[c].last_friction_rotational_impulse_A;
                rotational_impulses[b_id] -= contacts[c].last_friction_rotational_impulse_B;
                contacts[c].last_friction_rotational_impulse_A = rotational_impulse_friction_A;
                contacts[c].last_friction_rotational_impulse_B = rotational_impulse_friction_B;
                rotational_impulses[a_id] += contacts[c].last_friction_rotational_impulse_A;
                rotational_impulses[b_id] += contacts[c].last_friction_rotational_impulse_B;

                // ======== INTEGRATE ========

                if (!cluster.particles[a_id].is_static)
                {
                    Eigen::Vector3d impulse = impulses[a_id];
                    Eigen::Vector3d rotational_impulse = rotational_impulses[a_id];
                    Eigen::Vector3d impulse_offset = impulse_offsets[a_id];
                    Eigen::Vector3d rotational_impulse_offset = rotational_impulse_offsets[a_id];
                    FStates[a_id] = updateState(cluster.particles[a_id].future_state, t, impulse, impulse_offset, rotational_impulse, rotational_impulse_offset, &cluster.particles[a_id]);
                }
                if (!cluster.particles[b_id].is_static)
                {
                    Eigen::Vector3d impulse = impulses[b_id];
                    Eigen::Vector3d rotational_impulse = rotational_impulses[b_id];
                    Eigen::Vector3d impulse_offset = impulse_offsets[b_id];
                    Eigen::Vector3d rotational_impulse_offset = rotational_impulse_offsets[b_id];
                    FStates[b_id] = updateState(cluster.particles[b_id].future_state, t, impulse, impulse_offset, rotational_impulse, rotational_impulse_offset, &cluster.particles[a_id]);
                }
            }
            converged = true;
            for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                Eigen::Vector3d impulse_diff = impulses[p_i] - last_impulse[p_i];
                if (impulse_diff.norm() > 1e-4) {
                    converged = false;
                    break;
                }
                // printf("impulse: %f, %f, %f\n", impulses[p_i].x(), impulses[p_i].y(), impulses[p_i].z());
                // printf("impulse diff: %f, %f, %f\n", impulse_diff.x(), impulse_diff.y(), impulse_diff.z());
            }
            it++;
        }
        printf("Finish contact in %i iterations\n", it);
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        cluster.particles[p_i].future_state = FStates[p_i];
    }

    // common::Viewer view;
    // view.addParticleInterval(cluster.particles[0]);
    // view.addParticleInterval(cluster.particles[1]);
    // for (int c = 0; c < hits.size(); c++) {
    //     view.addEdge(common::Edge(hits[c].A, hits[c].B));
    // }
    // view.show();

    return true;
}
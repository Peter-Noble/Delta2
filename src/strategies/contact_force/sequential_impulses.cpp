#include "sequential_impulses.h"

#include "sequential_impulses_warm_start.h"
#include "LocalContact.h"
#include "update_state.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"
#include "../../model/friction.h"

#include "../../common/viewer.h"

#include <algorithm>
#include <atomic>

#include "tbb/tbb.h"

#include "../../globals.h"

// #include <ittnotify.h>

using namespace Delta2;
using namespace strategy;
using namespace Eigen;

class Spinlock
{
private:
    std::atomic_flag atomic_flag = ATOMIC_FLAG_INIT;

public:
    void lock()
    {
       while (atomic_flag.test_and_set(std::memory_order_acquire))
        {
        }
     
    }
    void unlock()
    {
        atomic_flag.clear(std::memory_order_release);
    }
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


State updateStateImmediate(State& future, double t, const Vector3d& impulse, const Vector3d& impulse_offset, const Vector3d& rotational_impulse, const Vector3d& rotational_impulse_offset, Particle* p) {
    State updated = future;

    updated.setVelocity(future.getVelocity() + impulse / p->getMass());
    updated.setTranslation(future.getTranslation() + impulse_offset / p->getMass() * t + impulse / p->getMass() * t);

    // https://link.springer.com/content/pdf/10.1007/s00707-013-0914-2.pdf

    Eigen::Vector3d angular_momentum = future.getAngularMomentum() + rotational_impulse;

    Eigen::Quaterniond rotation = future.getRotation();

    updated.setAngular(angular_momentum);

    Eigen::Matrix3d R = rotation.toRotationMatrix();
    Eigen::Matrix3d Iinv = R * p->getInverseInertiaMatrix() * R.transpose();
    Eigen::Vector3d omega = Iinv * (future.getAngularMomentum() + t * rotational_impulse_offset + t * rotational_impulse);

    Eigen::Quaterniond rot_delta = Delta2::common::exp(0.5 * omega * t);

    rotation = rot_delta * rotation;
    updated.setRotation(rotation);

    return updated;
}

SequentialImpulses::SequentialImpulses(FrictionStrategy& friction,
                                       common::Options& opt) :
                                       ContactForceStrategy(friction, opt) {

}

std::vector<LocalContact> create_local_contacts(std::vector<collision::Contact<double>>& hits) {
    std::vector<LocalContact> contacts;
    
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
        }
        if (!hits[c].p_b->is_static) {
            k_normal += im_B + rn_B.transpose() * ii_B * rn_B;
            k_tangent += im_B + rt_B.transpose() * ii_B * rt_B;
        }

        lc.mass_eff_normal = k_normal < 1e-6 ? 0.0 : 1.0 / k_normal;
        lc.mass_eff_tangent = k_tangent < 1e-6 ? 0.0 : 1.0 / k_tangent;

        lc.impulse = 0.0;
        lc.impulse_offset = 0.0;

        lc.last_friction_impulse = { 0.0, 0.0, 0.0 };
        lc.last_friction_rotational_impulse_A = { 0.0, 0.0, 0.0 };
        lc.last_friction_rotational_impulse_B = { 0.0, 0.0, 0.0 };

        contacts.push_back(lc);
    }
    return contacts;
}

void update_FStates(collision::Cluster& cluster,
                    std::vector<Eigen::Vector3d>& impulses,
                    std::vector<Eigen::Vector3d>& rotational_impulses,
                    std::vector<Eigen::Vector3d>& impulse_offsets,
                    std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                    std::vector<State>& FStates) {
    const double t = cluster.step_size;

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        if (!cluster.particles[p_i].is_static)
        {
            Eigen::Vector3d impulse = impulses[p_i];
            Eigen::Vector3d rotational_impulse = rotational_impulses[p_i];
            Eigen::Vector3d impulse_offset = impulse_offsets[p_i];
            Eigen::Vector3d rotational_impulse_offset = rotational_impulse_offsets[p_i];
            FStates[p_i] = updateState(cluster.particles[p_i].future_state, t, impulse, impulse_offset, rotational_impulse, rotational_impulse_offset, &cluster.particles[p_i]);
        }
    }
}

bool solve_contacts(collision::Cluster& cluster,
                    std::vector<collision::Contact<double>>& hits,
                    std::vector<LocalContact>& contacts,
                    std::vector<Eigen::Vector3d>& impulses,
                    std::vector<Eigen::Vector3d>& rotational_impulses,
                    std::vector<Eigen::Vector3d>& impulse_offsets,
                    std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                    std::vector<State>& FStates) {
    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    std::map<std::pair<int, int>, std::pair<int, Eigen::Vector3d>> average_contact_per_pair_first;

    const int max_iterations = globals::opt.sequential_impulse_total_iterations;
    const double t = cluster.step_size;

    bool all_d1_above = true;
    bool converged = false;
    int it = 0;
    while (!converged && it < max_iterations) {
        converged = true;
        all_d1_above = true;
        std::vector<Eigen::Vector3d> last_impulse;
        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            last_impulse.push_back(impulses[p_i]);
        }
        
        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

            uint32_t lower_id = std::min(a_id, b_id);
            uint32_t upper_id = std::max(a_id, b_id);

            Eigen::Matrix4d a_FState_transformation = FStates[a_id].getTransformation();
            Eigen::Matrix4d b_FState_transformation = FStates[b_id].getTransformation();
            Eigen::Matrix3d a_inverse_inertia_matrix = hits[c].p_a->getInverseInertiaMatrix();
            Eigen::Matrix3d b_inverse_inertia_matrix = hits[c].p_b->getInverseInertiaMatrix();
            Eigen::Vector3d a_FState_translation = FStates[a_id].getTranslation();
            Eigen::Vector3d b_FState_translation = FStates[b_id].getTranslation();

            Eigen::Vector3d n = contacts[c].global_normal.normalized();
            if (contacts[c].global_normal.norm() < 1e-12) {
                n = (b_FState_translation - a_FState_translation).normalized();
            }

            const double outer_search = hits[c].eps_a + hits[c].eps_b;
            const double outer = outer_search * 0.25;

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, a_FState_transformation);
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, b_FState_transformation);
            const double d1 = p1_B.dot(n) - p1_A.dot(n);

            all_d1_above &= d1 > 0.0;

            // double t0 = cluster.particles[a_id].current_state.getTime();
            // double t1 = FStates[a_id].getTime();
            // double ts = t1 - t0;

            // State FState_extrapolate_A = FStates[a_id].extrapolate(ts, hits[c].p_a->getInverseInertiaMatrix());
            // State FState_extrapolate_B = FStates[b_id].extrapolate(ts, hits[c].p_b->getInverseInertiaMatrix());

            // Eigen::Vector3d pf_A = common::transform(contacts[c].local_A, FState_extrapolate_A.getTransformation());
            // Eigen::Vector3d pf_B = common::transform(contacts[c].local_B, FState_extrapolate_B.getTransformation());

            // Eigen::Vector3d v1_A_ext = FStates[a_id].pointVelocity(p1_A, FState_extrapolate_A);
            // Eigen::Vector3d v1_B_ext = FStates[b_id].pointVelocity(p1_B, FState_extrapolate_B);

            Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, a_inverse_inertia_matrix);
            Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, b_inverse_inertia_matrix);
            
            const double v1 = v1_A.dot(n) - v1_B.dot(n);

            // Eigen::Vector3d v1_AB_ext = v1_A_ext - v1_B_ext;
            // const double v1 = v1_A_ext.dot(n) - v1_B_ext.dot(n);

            const double m = contacts[c].mass_eff_normal;

            // ======== CANCEL RELATIVE VELOCITY ALONG NORMAL ========

            const double delta_damp = common::lerp(0.5, 0.1, (double)it/(double)max_iterations);
            // v1 = im / m to cancel velocity
            const double delta_impulse_mag = d1 < outer ? v1 * m * delta_damp : -contacts[c].impulse;

            double last_impulse_mag = contacts[c].impulse;
            contacts[c].impulse = std::max(contacts[c].impulse + delta_impulse_mag, 0.0);

            Eigen::Vector3d total_impulse = n * contacts[c].impulse;
            // globals::logger.printf(3, "Contact %i Total impulse: %f, %f, %f\n", c, total_impulse.x(), total_impulse.y(), total_impulse.z());

            auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));
            if (it > 0) {
                if (average_contact_per_pair_first.find(key) == average_contact_per_pair_first.end()) {
                    // globals::logger.printf(3, "Searching pair %i, %i\n", key.first, key.second);
                    throw std::runtime_error("Key pair not found");
                }
                Eigen::Vector3d average = average_contact_per_pair_first[key].second / average_contact_per_pair_first[key].first;
                total_impulse = common::lerp(average, total_impulse, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
                total_impulse = d1 < outer ? n.dot(total_impulse) * n : Eigen::Vector3d({0, 0, 0});  // Keep direction of impulse
            }
            else {
                if (average_contact_per_pair_first.find(key) != average_contact_per_pair_first.end()) {
                    average_contact_per_pair_first[key].first++;
                    average_contact_per_pair_first[key].second += n * last_impulse_mag;
                }
                else {
                    average_contact_per_pair_first[key] = std::make_pair(1, n * last_impulse_mag);
                }
            }
            
            contacts[c].impulse = std::max(total_impulse.dot(n), 0.0);

            const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;

            Eigen::Vector3d update_impulse = n * update_impulse_mag;
            average_contact_per_pair_first[key].second += update_impulse;

            if (std::abs(update_impulse_mag) > 1e-8) {
                converged = false;
            }
            update_impulse *= -1;  //  The impulse is calculated A->B but applied B->A so flip it


            Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, a_FState_translation);
            Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, b_FState_translation);

            // ========  OFFSET ALONG NORMAL ========

            // const double lower_slop = std::max(outer * 0.25, outer - contacts[c].start_dist);
            // const double slop = common::lerp(lower_slop, outer * 0.25, 0.3);
            const double slop = outer * 0.5;
            const double penetration_depth = outer - d1;

            // d1 = i/m * t offset distance with impulse.  m * d1 / t = i
            double delta_impulse_offset_mag = d1 < outer ? delta_damp * m * std::max(0.0, penetration_depth - slop) / std::max(t, 1e-6) : -contacts[c].impulse_offset;

            if (std::isnan(delta_impulse_offset_mag)) {
                // throw std::runtime_error("offset mag is nan");
                delta_impulse_offset_mag = 0;
            }

            if (std::isinf(delta_impulse_offset_mag)) {
                // throw std::runtime_error("offset mag is inf");
                delta_impulse_offset_mag = 0;
            }

            double last_impulse_offset_mag = contacts[c].impulse_offset;
            contacts[c].impulse_offset = std::max(contacts[c].impulse_offset + delta_impulse_offset_mag, 0.0);
            const double update_impulse_offset_mag = contacts[c].impulse_offset - last_impulse_offset_mag;

            Eigen::Vector3d update_impulse_offset = -n * update_impulse_offset_mag;

            Eigen::Vector3d update_rotational_impulse_offset_A = model::calcTorque(update_impulse_offset, p1_A, a_FState_translation);
            Eigen::Vector3d update_rotational_impulse_offset_B = model::calcTorque(Eigen::Vector3d(-update_impulse_offset), p1_B, b_FState_translation);

            {
                // ======== INTEGRATE ========

                Eigen::Vector3d impulse_a;
                Eigen::Vector3d rotational_impulse_a;
                Eigen::Vector3d impulse_offset_a;
                Eigen::Vector3d rotational_impulse_offset_a;
                Eigen::Vector3d impulse_b;
                Eigen::Vector3d rotational_impulse_b;
                Eigen::Vector3d impulse_offset_b;
                Eigen::Vector3d rotational_impulse_offset_b;

                if (!cluster.particles[a_id].is_static)
                {
                    // #ifdef USETBB
                    // mxs[a_id]->lock();
                    // #endif

                    impulses[a_id] += update_impulse;
                    rotational_impulses[a_id] += update_rotational_impulse_A;
                    impulse_offsets[a_id] += update_impulse_offset;
                    rotational_impulse_offsets[a_id] += update_rotational_impulse_offset_A;

                    Eigen::Vector3d impulse_a = impulses[a_id];
                    Eigen::Vector3d rotational_impulse_a = rotational_impulses[a_id];
                    Eigen::Vector3d impulse_offset_a = impulse_offsets[a_id];
                    Eigen::Vector3d rotational_impulse_offset_a = rotational_impulse_offsets[a_id];

                    FStates[a_id] = updateState(cluster.particles[a_id].future_state, t, impulse_a, impulse_offset_a, rotational_impulse_a, rotational_impulse_offset_a, &cluster.particles[a_id]);
                    
                    if (FStates[a_id].getTransformation().hasNaN()) {
                        throw std::runtime_error("Nan");
                    }

                    // #ifdef USETBB
                    // mxs[a_id]->unlock();
                    // #endif
                }
                if (!cluster.particles[b_id].is_static)
                {
                    // #ifdef USETBB
                    // mxs[b_id]->lock();
                    // #endif

                    impulses[b_id] -= update_impulse;
                    rotational_impulses[b_id] += update_rotational_impulse_B;
                    impulse_offsets[b_id] -= update_impulse_offset;
                    rotational_impulse_offsets[b_id] += update_rotational_impulse_offset_B;

                    Eigen::Vector3d impulse_b = impulses[b_id];
                    Eigen::Vector3d rotational_impulse_b = rotational_impulses[b_id];
                    Eigen::Vector3d impulse_offset_b = impulse_offsets[b_id];
                    Eigen::Vector3d rotational_impulse_offset_b = rotational_impulse_offsets[b_id];

                    FStates[b_id] = updateState(cluster.particles[b_id].future_state, t, impulse_b, impulse_offset_b, rotational_impulse_b, rotational_impulse_offset_b, &cluster.particles[b_id]);
                    
                    if (FStates[b_id].getTransformation().hasNaN()) {
                        throw std::runtime_error("Nan");
                    }

                    // #ifdef USETBB
                    // mxs[b_id]->unlock();
                    // #endif
                }
            }
        }

        if (!all_d1_above && converged) {
            globals::logger.printf(4, "%i, Not all d1 above %i\n", cluster_id, it);
        }
        converged &= all_d1_above;

        it++;
    }

    if (!all_d1_above) {
        globals::logger.printf(3, "%i, returned due to not all d1 above\n", cluster_id);
        return false;
    }

    globals::logger.printf(3, "%i: Iterations used contact phase: %i\n", cluster_id, it);
    for (int c = 0; c < hits.size(); c++) {
        Eigen::Vector3d n = contacts[c].global_normal.normalized();
        Eigen::Vector3d total_impulse = n * contacts[c].impulse;
        globals::logger.printf(3, "Contact %i total_impulse: %f, %f, %f (%f)\n", c, total_impulse.x(), total_impulse.y(), total_impulse.z(), total_impulse.norm());
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        // assert(!impulses[p_i].hasNaN());
        // assert(impulses[p_i].norm() < 10000);
        // if (impulses[p_i].norm() > 10000) {
        //     globals::logger.printf(3, "Really big impulse\n");
        // }
        if (impulses[p_i].hasNaN()) {
            // // throw std::runtime_error("Nan!");
            globals::logger.printf(3, "Nan!\n");
            return false;
        }
    }
    return true;
}

bool solve_friction(collision::Cluster& cluster,
                    std::vector<collision::Contact<double>>& hits,
                    std::vector<LocalContact>& contacts,
                    std::vector<Eigen::Vector3d>& impulses,
                    std::vector<Eigen::Vector3d>& rotational_impulses,
                    std::vector<Eigen::Vector3d>& impulse_offsets,
                    std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                    std::vector<State>& FStates) {
    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }
    
    const int max_iterations = globals::opt.sequential_impulse_total_iterations;

    std::vector<Eigen::Vector3d> tangents;

    std::map<std::pair<int, int>, std::pair<int, Eigen::Vector3d>> average_contact_per_pair;
    
    std::vector<Eigen::Vector3d> updated_friction_total;
    for (int c = 0; c < hits.size(); c++) {
        updated_friction_total.push_back({0, 0, 0});
    }

    bool converged = false;
    int it = 0;
    while (!converged && it < max_iterations || it < 10) {
        converged = true;

        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

            Eigen::Vector3d n = contacts[c].global_normal.normalized();
            if (n.hasNaN()) {
                Eigen::Vector3d a_FState_translation = FStates[a_id].getTranslation();
                Eigen::Vector3d b_FState_translation = FStates[b_id].getTranslation();
                n = (b_FState_translation - a_FState_translation).normalized();
                globals::logger.printf(1, "Normal could not be computed in sequential impulse solver friction stage\n");
            }

            double t0 = cluster.particles[a_id].current_state.getTime();
            double t1 = FStates[a_id].getTime();
            double ts = t1 - t0;

            State FState_extrapolate_A = FStates[a_id].extrapolate(ts, hits[c].p_a->getInverseInertiaMatrix());
            State FState_extrapolate_B = FStates[b_id].extrapolate(ts, hits[c].p_b->getInverseInertiaMatrix());

            Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
            Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

            Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());
            Eigen::Vector3d v1_AB = v1_A - v1_B;
            Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
            
            // Could alternatively be used to compute the velocity at this point as the finite difference between this point now and at the last time step
            // Eigen::Vector3d v1_A_fd = (p1_A - p0_A) / ts;
            // Eigen::Vector3d v1_B_fd = (p1_B - p0_B) / ts;
            // Eigen::Vector3d v1_AB_fd = v1_A_fd - v1_B_fd;
            // Eigen::Vector3d v_tangent_fd = v1_AB_fd - v1_AB_fd.dot(n) * n;

            // A third alternative.  Project the state out into the future and take a finite difference with that.
            // Eigen::Vector3d v1_A_ext = FStates[a_id].pointVelocity(p1_A, FState_extrapolate_A);
            // Eigen::Vector3d v1_B_ext = FStates[b_id].pointVelocity(p1_B, FState_extrapolate_B);
            // Eigen::Vector3d v1_AB_ext = v1_A_ext - v1_B_ext;
            // Eigen::Vector3d v_tangent_ext = v1_AB_ext - v1_AB_ext.dot(n) * n;

            Eigen::Vector3d v_tangent = v_tangent_pt;
            Eigen::Vector3d tangent = v_tangent.normalized();

            if (it == 0) {
                globals::logger.printf(3, "%i: %i v_tangent: %f, %f, %f (normal %f, tangent %f)\n", cluster_id, c, v_tangent.x(), v_tangent.y(), v_tangent.z(), contacts[c].impulse, contacts[c].last_friction_impulse.norm());
                globals::logger.printf(3, "%i: %i v1_AB: %f, %f, %f first friction\n", cluster_id, c, v1_AB.x(), v1_AB.y(), v1_AB.z());
                globals::logger.printf(3, "%i: %i n: %f, %f, %f first friction\n", cluster_id, c, n.x(), n.y(), n.z());
            }
            
            double mu = std::min(cluster.particles[a_id].friction_coeff, cluster.particles[b_id].friction_coeff);
            double max_impulse = (contacts[c].impulse + contacts[c].impulse_offset) * mu;

            Eigen::Vector3d r_A = hits[c].A - hits[c].p_a->future_state.getTranslation();
            Eigen::Vector3d r_B = hits[c].B - hits[c].p_b->future_state.getTranslation();

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

            if (k_tangent < 1e-6 && max_impulse > 0) {
                globals::logger.printf(1, "Extremely small k_tangent for non-zero impulse");
            }
            double m_inv_effective_tangent = k_tangent < 1e-6 ? 0.0 : 1.0 / k_tangent;

            const double friction_damp = common::lerp(0.1, 0.01, (double)it/(double)max_iterations);
            Eigen::Vector3d delta_friction_impulse = -v_tangent * m_inv_effective_tangent * friction_damp;
            // Eigen::Vector3d delta_friction_impulse_fd = -v_tangent_fd * m_inv_effective_tangent * friction_damp;

            Eigen::Vector3d friction_impulse_total = contacts[c].last_friction_impulse + delta_friction_impulse;
            Eigen::Vector3d friction_impulse_total_pre_avg = friction_impulse_total;

            if (it > 0) {
                auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));
                Eigen::Vector3d average = average_contact_per_pair[key].second / average_contact_per_pair[key].first;
                friction_impulse_total = common::lerp(average, friction_impulse_total, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
                friction_impulse_total = -tangent * std::max(0.0, (-tangent).dot(friction_impulse_total));
            }
            friction_impulse_total = friction_impulse_total_pre_avg.normalized() * std::clamp(friction_impulse_total.norm(), 0.0, max_impulse);

            updated_friction_total[c] = friction_impulse_total;

            if ((contacts[c].last_friction_impulse - updated_friction_total[c]).norm() > 1e-8) {
                converged = false;
            }

            if (updated_friction_total[c].hasNaN()) {
                globals::logger.printf(1, "Updated friction impulse has nan\n");
                globals::logger.printf(1, "c: %i\n", c);
                globals::logger.printf(1, FStates[a_id].getTransformation().hasNaN() ? "FStates a_id has nan" : "FStates a_id is valid\n");
                globals::logger.printf(1, FStates[b_id].getTransformation().hasNaN() ? "FStates b_id has nan" : "FStates b_id is valid\n");
                globals::logger.printf(1, "contacts[c].local_A: %f, %f, %f\n", contacts[c].local_A.x(), contacts[c].local_A.y(), contacts[c].local_A.z());
                globals::logger.printf(1, "contacts[c].local_B: %f, %f, %f\n", contacts[c].local_B.x(), contacts[c].local_B.y(), contacts[c].local_B.z());
                globals::logger.printf(1, "p1_A: %f, %f, %f\n", p1_A.x(), p1_A.y(), p1_A.z());
                globals::logger.printf(1, "p1_B: %f, %f, %f\n", p1_B.x(), p1_B.y(), p1_B.z());
                globals::logger.printf(1, "v1_AB: %f, %f, %f\n", v1_AB.x(), v1_AB.y(), v1_AB.z());
                globals::logger.printf(1, "rt_A: %f, %f, %f\n", rt_A.x(), rt_A.y(), rt_A.z());
                globals::logger.printf(1, "rt_B: %f, %f, %f\n", rt_B.x(), rt_B.y(), rt_B.z());
                globals::logger.printf(1, "tangent: %f, %f, %f\n", tangent.x(), tangent.y(), tangent.z());
                globals::logger.printf(1, "k_tangent: %f\n", k_tangent);
                globals::logger.printf(1, "m_inv_effective_tangent: %f\n", m_inv_effective_tangent);
                globals::logger.printf(1, "delta_friction_impulse: %f, %f, %f\n", delta_friction_impulse.x(), delta_friction_impulse.y(), delta_friction_impulse.z());
                // throw std::runtime_error("Nans");
            }
        }

        average_contact_per_pair.clear();

        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
            auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));

            if (average_contact_per_pair.find(key) != average_contact_per_pair.end()) {
                average_contact_per_pair[key].first++;
                average_contact_per_pair[key].second += updated_friction_total[c];
            }
            else {
                average_contact_per_pair[key] = std::make_pair(1, updated_friction_total[c]);
            }

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

            Eigen::Vector3d update_friction_impulse = (updated_friction_total[c] - contacts[c].last_friction_impulse) / average_contact_per_pair[key].first;
            contacts[c].last_friction_impulse = updated_friction_total[c];
            impulses[a_id] += update_friction_impulse;
            impulses[b_id] -= update_friction_impulse;

            // impulse_offsets[a_id] += update_friction_impulse;
            // impulse_offsets[b_id] -= update_friction_impulse;

            Eigen::Vector3d rotational_impulse_friction_A = model::calcTorque(contacts[c].last_friction_impulse, p1_A, FStates[a_id].getTranslation());
            Eigen::Vector3d rotational_impulse_friction_B = model::calcTorque(Eigen::Vector3d(-contacts[c].last_friction_impulse), p1_B, FStates[b_id].getTranslation());

            rotational_impulses[a_id] -= contacts[c].last_friction_rotational_impulse_A;
            rotational_impulses[b_id] -= contacts[c].last_friction_rotational_impulse_B;
            // rotational_impulse_offsets[a_id] -= contacts[c].last_friction_rotational_impulse_A;
            // rotational_impulse_offsets[b_id] -= contacts[c].last_friction_rotational_impulse_B;
            contacts[c].last_friction_rotational_impulse_A = rotational_impulse_friction_A;
            contacts[c].last_friction_rotational_impulse_B = rotational_impulse_friction_B;
            rotational_impulses[a_id] += contacts[c].last_friction_rotational_impulse_A;
            rotational_impulses[b_id] += contacts[c].last_friction_rotational_impulse_B;
            // rotational_impulse_offsets[a_id] += contacts[c].last_friction_rotational_impulse_A;
            // rotational_impulse_offsets[b_id] += contacts[c].last_friction_rotational_impulse_B;
        }

        // ======== INTEGRATE ========
        update_FStates(cluster, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates);
        it++;
    }

    // if (it == max_iterations) {
    //     return false;  // TODO this used to cause some very odd errors in the tiny_tower scenario. Still the case?
    // }

    for (int c = 0; c < hits.size(); c++) {
        uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
        uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

        Eigen::Vector3d n = contacts[c].global_normal.normalized();

        Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
        Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());
        
        Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
        Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

        double t0 = cluster.particles[a_id].current_state.getTime();
        double t1 = FStates[a_id].getTime();
        double ts = t1 - t0;

        Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
        Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());
        Eigen::Vector3d v1_AB = v1_A - v1_B;
        
        Eigen::Vector3d v1_A_fd = (p1_A - p0_A) / ts;
        Eigen::Vector3d v1_B_fd = (p1_B - p0_B) / ts;
        Eigen::Vector3d v1_AB_fd = v1_A_fd - v1_B_fd;

        Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
        Eigen::Vector3d v_tangent_fd = v1_AB_fd - v1_AB_fd.dot(n) * n;

        Eigen::Vector3d v_tangent = v_tangent_pt;
        Eigen::Vector3d tangent = v_tangent.normalized();

        globals::logger.printf(4, "%i: %i v_tangent: %f, %f, %f, (normal %f, tangent %f)\n", cluster_id, c, v_tangent.x(), v_tangent.y(), v_tangent.z(), contacts[c].impulse, contacts[c].last_friction_impulse.norm());
        globals::logger.printf(4, "%i: %i v1_AB: %f, %f, %f last friction\n", cluster_id, c, v1_AB.x(), v1_AB.y(), v1_AB.z());
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        assert(!impulses[p_i].hasNaN());
        if (impulses[p_i].hasNaN()) {
            globals::logger.printf(3, "Nan!\n");
            return false;
        }
    }
    
    globals::logger.printf(3, "%i: Iterations used friction phase: %i\n", cluster_id, it);

    for (int c = 0; c < hits.size(); c++) {
        globals::logger.printf(3, "Contact %i last_impulse_friction: %f, %f, %f (%f)\n", c, contacts[c].last_friction_impulse.x(), contacts[c].last_friction_impulse.y(), contacts[c].last_friction_impulse.z(), contacts[c].last_friction_impulse.norm());
    }

    return true;
}

bool SequentialImpulses::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits, bool allow_fail) {
    std::sort(hits.begin(), hits.end());

    std::vector<LocalContact> contacts = create_local_contacts(hits);

    std::vector<Eigen::Vector3d> impulses;
    std::vector<Eigen::Vector3d> rotational_impulses;
    std::vector<Eigen::Vector3d> impulse_offsets;
    std::vector<Eigen::Vector3d> rotational_impulse_offsets;

    std::vector<State> FStates;

    const double t = cluster.step_size;

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        impulses.push_back(external_force(cluster.particles[p_i]) * t);
        rotational_impulses.push_back({0, 0, 0});
        impulse_offsets.push_back({0, 0, 0});
        rotational_impulse_offsets.push_back({0, 0, 0});
        FStates.push_back(cluster.particles[p_i].future_state);
    }

    // _warm_start.apply(contacts, hits);

    // Apply impulses from warm start
    for (int c = 0; c < hits.size(); c++) {
        uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
        uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

        Eigen::Vector3d n = contacts[c].global_normal.normalized();
        if (n.hasNaN()) {
            continue;
        }

        Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
        Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
        const double d1 = p1_B.dot(n) - p1_A.dot(n);

        Eigen::Vector3d update_impulse = -n * contacts[c].impulse;

        impulses[a_id] += update_impulse;
        impulses[b_id] -= update_impulse;

        Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, FStates[a_id].getTranslation());
        Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, FStates[b_id].getTranslation());

        rotational_impulses[a_id] += update_rotational_impulse_A;
        rotational_impulses[b_id] += update_rotational_impulse_B;
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        if (!cluster.particles[p_i].is_static) {
            Eigen::Vector3d velocity = FStates[p_i].getVelocity();
            globals::logger.printf(3, "Pre solve velocity for %i: %f, %f, %f\n", p_i, velocity.x(), velocity.y(), velocity.z());
            globals::logger.printf(3, "impulses[p_i] %f, %f %f\n", impulses[p_i].x(), impulses[p_i].y(), impulses[p_i].z());
        }
    }

    // Make a first prediction at the future state
    update_FStates(cluster, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates);


    if (hits.size() > 0) {
        if (!solve_contacts(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates)) return false;

        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            if (!cluster.particles[p_i].is_static) {
                Eigen::Vector3d velocity = FStates[p_i].getVelocity();
                globals::logger.printf(3, "Post solve velocity for %i: %f, %f, %f\n", p_i, velocity.x(), velocity.y(), velocity.z());
            }
        }   

        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
            
            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

            Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

            Eigen::Vector3d v1_AB = v1_A - v1_B;

            globals::logger.printf(3, "%i v1_AB: %f, %f, %f post 1st contact solve\n", c, v1_AB.x(), v1_AB.y(), v1_AB.z());
        }

        for (int inner = 0; inner < globals::opt.sequential_impulse_inner_iterations; inner++) {
            if (!solve_friction(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates)) return false;

            for (int p = 0; p < cluster.particles.size(); p++) {
                if (!cluster.particles[p].is_static) {
                    Eigen::Vector3d impulse = impulses[p];
                    Eigen::Vector3d rotational_impulse = rotational_impulses[p];
                    globals::logger.printf(3, "pf Impulse: %f, %f, %f\n", impulse.x(), impulse.y(), impulse.z());
                    globals::logger.printf(3, "pf rotational_impulse: %f, %f, %f\n", rotational_impulse.x(), rotational_impulse.y(), rotational_impulse.z());
                    globals::logger.printf(3, "pf Velocity %f, %f, %f\n", FStates[p].getVelocity().x(), FStates[p].getVelocity().y(), FStates[p].getVelocity().z());
                }
            }

            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
                Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());
                
                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

                double t0 = cluster.particles[a_id].current_state.getTime();
                double t1 = FStates[a_id].getTime();
                double ts = t1 - t0;
                
                Eigen::Vector3d v1_A_fd = (p1_A - p0_A) / ts;
                Eigen::Vector3d v1_B_fd = (p1_B - p0_B) / ts;

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

                // globals::logger.printf(4, "ts: %f\n", ts);
                // globals::logger.printf(4, "v1_A_fd: (%f, %f, %f) v1_A: (%f, %f, %f)\n", v1_A_fd.x(), v1_A_fd.y(), v1_A_fd.z(), v1_A.x(), v1_A.y(), v1_A.z());
                // globals::logger.printf(4, "v1_B_fd: (%f, %f, %f) v1_B: (%f, %f, %f)\n", v1_B_fd.x(), v1_B_fd.y(), v1_B_fd.z(), v1_B.x(), v1_B.y(), v1_B.z());

                Eigen::Vector3d v1_AB = v1_A - v1_B;
                Eigen::Vector3d v1_AB_fd = v1_A_fd - v1_B_fd;
                
                Eigen::Vector3d n = contacts[c].global_normal.normalized();

                Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
                Eigen::Vector3d v_tangent_fd = v1_AB_fd - v1_AB_fd.dot(n) * n;

                Eigen::Vector3d v_tangent = v_tangent_pt;
                Eigen::Vector3d tangent = v_tangent.normalized();

                globals::logger.printf(3, "%i v_tangent: %f, %f, %f, (normal %f, tangent %f) post friction\n", c, v_tangent.x(), v_tangent.y(), v_tangent.z(), contacts[c].impulse, contacts[c].last_friction_impulse.norm());
            }

            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
                
                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

                Eigen::Vector3d v1_AB = v1_A - v1_B;

                globals::logger.printf(3, "%i v1_AB: %f, %f, %f post friction solve\n", c, v1_AB.x(), v1_AB.y(), v1_AB.z());
            }

            if (!solve_contacts(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates)) return false;
            
            for (int p = 0; p < cluster.particles.size(); p++) {
                if (!cluster.particles[p].is_static) {
                    Eigen::Vector3d impulse = impulses[p];
                    Eigen::Vector3d rotational_impulse = rotational_impulses[p];
                    globals::logger.printf(3, "p2 Impulse: %f, %f, %f\n", impulse.x(), impulse.y(), impulse.z());
                    globals::logger.printf(3, "p2 rotational_impulse: %f, %f, %f\n", rotational_impulse.x(), rotational_impulse.y(), rotational_impulse.z());
                    globals::logger.printf(3, "p2 Velocity %f, %f, %f\n", FStates[p].getVelocity().x(), FStates[p].getVelocity().y(), FStates[p].getVelocity().z());
                }
            }

            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

                Eigen::Vector3d v1_AB = v1_A - v1_B;

                globals::logger.printf(3, "%i v1_AB: %f, %f, %f post 2nd contact solve\n", c, v1_AB.x(), v1_AB.y(), v1_AB.z());
            }
        }

        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
                globals::contact_draws.push_back(std::make_tuple(contacts[c].global_centre, contacts[c].global_normal.normalized() * contacts[c].impulse, contacts[c].last_friction_impulse));
            }
        }
    }

    // Apply integration
    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        cluster.particles[p_i].future_state = FStates[p_i];
    }

    for (int c = 0; c < hits.size(); c++) {
        uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
        uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

        Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
        Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());
        
        Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, cluster.particles[a_id].future_state.getTransformation());
        Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, cluster.particles[b_id].future_state.getTransformation());

        double t0 = cluster.particles[a_id].current_state.getTime();
        double t1 = cluster.particles[a_id].future_state.getTime();
        double ts = t1 - t0;
        
        Eigen::Vector3d v1_A_fd = (p1_A - p0_A) / ts;
        Eigen::Vector3d v1_B_fd = (p1_B - p0_B) / ts;

        Eigen::Vector3d v1_A = cluster.particles[a_id].future_state.pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
        Eigen::Vector3d v1_B = cluster.particles[b_id].future_state.pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

        // globals::logger.printf(4, "ts: %f\n", ts);
        // globals::logger.printf(4, "v1_A_fd: (%f, %f, %f) v1_A: (%f, %f, %f)\n", v1_A_fd.x(), v1_A_fd.y(), v1_A_fd.z(), v1_A.x(), v1_A.y(), v1_A.z());
        // globals::logger.printf(4, "v1_B_fd: (%f, %f, %f) v1_B: (%f, %f, %f)\n", v1_B_fd.x(), v1_B_fd.y(), v1_B_fd.z(), v1_B.x(), v1_B.y(), v1_B.z());

        Eigen::Vector3d v1_AB = v1_A - v1_B;
        Eigen::Vector3d v1_AB_fd = v1_A_fd - v1_B_fd;
           
        Eigen::Vector3d n = contacts[c].global_normal.normalized();

        Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
        Eigen::Vector3d v_tangent_fd = v1_AB_fd - v1_AB_fd.dot(n) * n;

        Eigen::Vector3d v_tangent = v_tangent_pt;
        Eigen::Vector3d tangent = v_tangent.normalized();

        globals::logger.printf(3, "%i v_tangent: %f, %f, %f, (normal %f, tangent %f) post integrate\n", c, v_tangent.x(), v_tangent.y(), v_tangent.z(), contacts[c].impulse, contacts[c].last_friction_impulse.norm());
    }

    // _warm_start.update(contacts, hits);

    return true;
}

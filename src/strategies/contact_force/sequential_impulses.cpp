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

        lc.last_friction_impulse_fd = { 0.0, 0.0, 0.0 };
        lc.last_friction_rotational_impulse_A_fd = { 0.0, 0.0, 0.0 };
        lc.last_friction_rotational_impulse_B_fd = { 0.0, 0.0, 0.0 };

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

    std::map<std::pair<int, int>, std::pair<int, Eigen::Vector3d>> contacts_per_pair_first;

    const int max_iterations = globals::opt.sequential_impulse_total_iterations;
    const double t = cluster.step_size;

    bool all_d1_above = true;
    bool converged = false;
    int it = 0;
    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_impulses_task);
    while (!converged && it < max_iterations) {
        converged = true;
        all_d1_above = true;
        if (it < max_iterations - 1) {
            // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_iteration_task);
        }
        else {
            // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_last_iteration_task);
        }
        std::vector<Eigen::Vector3d> last_impulse;
        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            last_impulse.push_back(impulses[p_i]);
        }
        
        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

            uint32_t lower_id = std::min(a_id, b_id);
            uint32_t upper_id = std::max(a_id, b_id);

            Eigen::Matrix4d a_FState_transformation;
            Eigen::Matrix4d b_FState_transformation;
            Eigen::Matrix3d a_inverse_inertia_matrix;
            Eigen::Matrix3d b_inverse_inertia_matrix;
            Eigen::Vector3d a_FState_translation;
            Eigen::Vector3d b_FState_translation;
            
            {
                // #ifdef USETBB
                // // particle_access_lock.lock();
                // mxs[lower_id]->lock();
                // mxs[upper_id]->lock();
                // #endif

                a_FState_transformation = FStates[a_id].getTransformation();
                b_FState_transformation = FStates[b_id].getTransformation();
                a_inverse_inertia_matrix = hits[c].p_a->getInverseInertiaMatrix();
                b_inverse_inertia_matrix = hits[c].p_b->getInverseInertiaMatrix();
                a_FState_translation = FStates[a_id].getTranslation();
                b_FState_translation = FStates[b_id].getTranslation();

                // #ifdef USETBB
                // mxs[lower_id]->unlock();
                // mxs[upper_id]->unlock();
                // // particle_access_lock.unlock();
                // #endif
            }

            Eigen::Vector3d n = contacts[c].global_normal.normalized();
            if (n.hasNaN()) {
                n = (b_FState_translation - a_FState_translation).normalized();
            }

            const double outer_search = hits[c].eps_a + hits[c].eps_b;
            const double outer = outer_search * 0.25;

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, a_FState_transformation);
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, b_FState_transformation);
            const double d1 = p1_B.dot(n) - p1_A.dot(n);

            all_d1_above &= d1 > 0.0;

            Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, a_inverse_inertia_matrix);
            Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, b_inverse_inertia_matrix);
            const double v1 = v1_A.dot(n) - v1_B.dot(n);

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
                if (contacts_per_pair_first.find(key) == contacts_per_pair_first.end()) {
                    // globals::logger.printf(3, "Searching pair %i, %i\n", key.first, key.second);
                    throw std::runtime_error("Key pair not found");
                }
                Eigen::Vector3d average = contacts_per_pair_first[key].second / contacts_per_pair_first[key].first;
                // globals::logger.printf(3, "Contact %i Average: %f, %f, %f\n", c, average.x(), average.y(), average.z());
                total_impulse = common::lerp(average, total_impulse, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
                total_impulse = n.dot(total_impulse) * n;  // Keep direction of impulse
                // globals::logger.printf(3, "Contact %i total_impulse: %f, %f, %f, lerp: %f\n", c, total_impulse.x(), total_impulse.y(), total_impulse.z(), std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
            }
            else {
                if (contacts_per_pair_first.find(key) != contacts_per_pair_first.end()) {
                    contacts_per_pair_first[key].first++;
                    contacts_per_pair_first[key].second += n * last_impulse_mag;
                }
                else {
                    // globals::logger.printf(3, "Adding new pair %i, %i\n", key.first, key.second);
                    contacts_per_pair_first[key] = std::make_pair(1, n * last_impulse_mag);
                }
                // globals::logger.printf(3, "per pair total: %f, %f, %f\n", contacts_per_pair_first[key].second.x(), contacts_per_pair_first[key].second.y(), contacts_per_pair_first[key].second.z());
            }
            
            contacts[c].impulse = std::max(total_impulse.dot(n), 0.0);

            const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;
            // const double update_impulse_mag = total_impulse.dot(n) - last_impulse_mag;

            Eigen::Vector3d update_impulse = n * update_impulse_mag;
            contacts_per_pair_first[key].second += update_impulse;

            // globals::logger.printf(3, "update_impulse_mag * 1000: %f\n", update_impulse_mag * 1000);
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
                        // throw std::runtime_error("Nan");
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
                        // throw std::runtime_error("Nan");
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

    double unsolved_velocity_total = 0.0;
    double last_velocity_total = 0.0;

    std::map<std::pair<int, int>, std::pair<int, Eigen::Vector3d>> contacts_per_pair;
    
    std::vector<Eigen::Vector3d> updated_friction_total;
    std::vector<Eigen::Vector3d> updated_friction_total_fd;
    for (int c = 0; c < hits.size(); c++) {
        updated_friction_total.push_back({0, 0, 0});
        updated_friction_total_fd.push_back({0, 0, 0});
    }

    bool converged = false;
    int it = 0;
    while (!converged && it < max_iterations || it < 10) {
        last_velocity_total = 0.0;
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

            Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
            Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

            Eigen::Vector3d v1_A_fd = (p1_A - p0_A) / ts;
            Eigen::Vector3d v1_B_fd = (p1_B - p0_B) / ts;

            Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
            Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

            // globals::logger.printf(4, "ts: %f\n", ts);
            // globals::logger.printf(4, "v1_A_fd: (%f, %f, %f) v1_A: (%f, %f, %f)\n", v1_A_fd.x(), v1_A_fd.y(), v1_A_fd.z(), v1_A.x(), v1_A.y(), v1_A.z());
            // globals::logger.printf(4, "v1_B_fd: (%f, %f, %f) v1_B: (%f, %f, %f)\n", v1_B_fd.x(), v1_B_fd.y(), v1_B_fd.z(), v1_B.x(), v1_B.y(), v1_B.z());

            Eigen::Vector3d v1_AB = v1_A - v1_B;
            Eigen::Vector3d v1_AB_fd = v1_A_fd - v1_B_fd;

            Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
            Eigen::Vector3d v_tangent_fd = v1_AB_fd - v1_AB_fd.dot(n) * n;

            Eigen::Vector3d v_tangent = v_tangent_pt;
            Eigen::Vector3d tangent = v_tangent.normalized();

            last_velocity_total += v_tangent_fd.norm();
            if (it == 0) {
                unsolved_velocity_total += v_tangent_fd.norm();
            //     globals::logger.printf(2, "v_tangent_fd norm: (%f), %f, %f, %f\n", v_tangent_fd.norm(), v_tangent_fd.x(), v_tangent_fd.y(), v_tangent_fd.z());
            //     Eigen::Vector3d tangent = v_tangent_fd.normalized();
            //     tangents.push_back(tangent);
            }
            // Eigen::Vector3d tangent = tangents[c];

            // v_tangent = v_tangent.dot(tangents[c]) * tangents[c];
            // v_tangent_fd = v_tangent_fd.dot(tangents[c]) * tangents[c];

            // globals::logger.printf(4, "v_tangent: %f, %f, %f, v_tangent_fd: %f, %f, %f\n", v_tangent.x(), v_tangent.y(), v_tangent.z(), v_tangent_fd.x(), v_tangent_fd.y(), v_tangent_fd.z());

            double mu = std::min(cluster.particles[a_id].friction_coeff, cluster.particles[b_id].friction_coeff);
            double max_impulse = contacts[c].impulse * mu;

            // globals::logger.printf(3, "it: %i, c: %i, max_impulse: %.4f, v_tangent: %.4f\n", it, c, max_impulse, v_tangent.norm());
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
            Eigen::Vector3d delta_friction_impulse_fd = -v_tangent_fd * m_inv_effective_tangent * friction_damp;

            // globals::logger.printf(4, "max_impulse: %f, delta_friction_impulse: %.4f, %.4f, %.4f\n", max_impulse, delta_friction_impulse.x(), delta_friction_impulse.y(), delta_friction_impulse.z());

            // double zero_bias = common::lerp(max_impulse * 0.1, 0.0, std::min(1.0, (double)it/10.0));
            double zero_bias = 0.0;
            Eigen::Vector3d friction_impulse_total = contacts[c].last_friction_impulse + delta_friction_impulse;
            Eigen::Vector3d friction_impulse_total_pre_avg = friction_impulse_total;

            if (it > 0) {
                auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));
                Eigen::Vector3d average = contacts_per_pair[key].second / contacts_per_pair[key].first;
                friction_impulse_total = common::lerp(average, friction_impulse_total, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
                friction_impulse_total = -tangent * std::max(0.0, (-tangent).dot(friction_impulse_total));
            }
            friction_impulse_total = friction_impulse_total_pre_avg.normalized() * std::clamp(friction_impulse_total.norm() - zero_bias, 0.0, max_impulse);

            updated_friction_total[c] = friction_impulse_total;

            Eigen::Vector3d friction_impulse_total_fd = contacts[c].last_friction_impulse_fd + delta_friction_impulse_fd;
            friction_impulse_total_fd = friction_impulse_total_fd.normalized() * std::clamp(friction_impulse_total_fd.norm(), 0.0, max_impulse);

            updated_friction_total_fd[c] = friction_impulse_total_fd;

            // globals::logger.printf(4, "lfi: %f, %f, %f, uft: %f, %f, %f\n", contacts[c].last_friction_impulse.x(), contacts[c].last_friction_impulse.y(), contacts[c].last_friction_impulse.z(), updated_friction_total[c].x(), updated_friction_total[c].y(), updated_friction_total[c].z());
            if ((contacts[c].last_friction_impulse - updated_friction_total[c]).norm() > 1e-8) {  // TODO This was originally 1e-8
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

        contacts_per_pair.clear();

        for (int c = 0; c < hits.size(); c++) {
            uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
            uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
            auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));

            if (contacts_per_pair.find(key) != contacts_per_pair.end()) {
                // globals::logger.printf(2, "Adding new pair %i, %i\n", key.first, key.second);
                contacts_per_pair[key].first++;
                contacts_per_pair[key].second += updated_friction_total[c];
            }
            else {
                contacts_per_pair[key] = std::make_pair(1, updated_friction_total[c]);
            }

            Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
            Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

            Eigen::Vector3d update_friction_impulse = (updated_friction_total[c] - contacts[c].last_friction_impulse) / contacts_per_pair[key].first;
            contacts[c].last_friction_impulse = updated_friction_total[c];
            impulses[a_id] += update_friction_impulse;
            impulses[b_id] -= update_friction_impulse;

            Eigen::Vector3d update_friction_impulse_fd = (updated_friction_total_fd[c] - contacts[c].last_friction_impulse_fd) / contacts_per_pair[key].first;
            contacts[c].last_friction_impulse_fd = updated_friction_total_fd[c];
            impulse_offsets[a_id] += update_friction_impulse_fd;
            impulse_offsets[b_id] -= update_friction_impulse_fd;

            Eigen::Vector3d rotational_impulse_friction_A = model::calcTorque(contacts[c].last_friction_impulse, p1_A, FStates[a_id].getTranslation());
            Eigen::Vector3d rotational_impulse_friction_B = model::calcTorque(Eigen::Vector3d(-contacts[c].last_friction_impulse), p1_B, FStates[b_id].getTranslation());

            rotational_impulses[a_id] -= contacts[c].last_friction_rotational_impulse_A;
            rotational_impulses[b_id] -= contacts[c].last_friction_rotational_impulse_B;
            contacts[c].last_friction_rotational_impulse_A = rotational_impulse_friction_A;
            contacts[c].last_friction_rotational_impulse_B = rotational_impulse_friction_B;
            rotational_impulses[a_id] += contacts[c].last_friction_rotational_impulse_A;
            rotational_impulses[b_id] += contacts[c].last_friction_rotational_impulse_B;

            Eigen::Vector3d rotational_impulse_friction_A_fd = model::calcTorque(contacts[c].last_friction_impulse_fd, p1_A, FStates[a_id].getTranslation());
            Eigen::Vector3d rotational_impulse_friction_B_fd = model::calcTorque(Eigen::Vector3d(-contacts[c].last_friction_impulse_fd), p1_B, FStates[b_id].getTranslation());

            rotational_impulse_offsets[a_id] -= contacts[c].last_friction_rotational_impulse_A_fd;
            rotational_impulse_offsets[b_id] -= contacts[c].last_friction_rotational_impulse_B_fd;
            contacts[c].last_friction_rotational_impulse_A_fd = rotational_impulse_friction_A_fd;
            contacts[c].last_friction_rotational_impulse_B_fd = rotational_impulse_friction_B_fd;
            rotational_impulse_offsets[a_id] += contacts[c].last_friction_rotational_impulse_A_fd;
            rotational_impulse_offsets[b_id] += contacts[c].last_friction_rotational_impulse_B_fd;
        }

        // ======== INTEGRATE ========
        update_FStates(cluster, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates);
        it++;
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        assert(!impulses[p_i].hasNaN());
        // assert(impulses[p_i].norm() < 10000);
        // if (impulses[p_i].norm() > 10000) {
        //     globals::logger.printf(3, "Really big impulse\n");
        // }
        if (impulses[p_i].hasNaN()) {
            // throw std::runtime_error("Nan!");
            globals::logger.printf(3, "Nan!\n");
            return false;
        }
    }
    
    globals::logger.printf(3, "%i: Iterations used friction phase: %i\n", cluster_id, it);
    globals::logger.printf(3, "Initial velocity lose: %f, end velocity lose: %f\n", unsolved_velocity_total, last_velocity_total);

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

    _warm_start.update(contacts, hits);

    return true;
}

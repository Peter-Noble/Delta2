#include "sequential_impulses.h"

#include "sequential_impulses_warm_start.h"
#include "LocalContact.h"
#include "update_state.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"
#include "../../model/friction.h"

#include "../../common/viewer.h"

#include "../../collision_detection/colour_hits.h"

#include <algorithm>
#include <atomic>
#include <chrono>

#include "tbb/tbb.h"

#include "../../globals.h"

// #include <ittnotify.h>

using namespace Delta2;
using namespace strategy;
using namespace Eigen;

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::microseconds;
using namespace std::literals::chrono_literals;
using std::this_thread::sleep_for;

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

void contact_iteration(collision::ContactBundle& bundle,
                       int it,
                       collision::Cluster& cluster,
                       std::vector<collision::Contact<double>>& hits,
                       std::vector<LocalContact>& contacts,
                       std::vector<Eigen::Vector3d>& impulses,
                       std::vector<Eigen::Vector3d>& rotational_impulses,
                       std::vector<Eigen::Vector3d>& impulse_offsets,
                       std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                       std::vector<State>& FStates,
                       std::atomic_bool& all_d1_above,
                       std::atomic_bool& converged) {
    const int max_iterations = globals::opt.sequential_impulse_total_iterations;
    const double t = cluster.step_size;
                        
    const uint32_t a_id = bundle.lower;
    const uint32_t b_id = bundle.upper;

    const Eigen::Matrix3d a_inverse_inertia_matrix = cluster.particles[a_id].getInverseInertiaMatrix();
    const Eigen::Matrix3d b_inverse_inertia_matrix = cluster.particles[b_id].getInverseInertiaMatrix();

    for (int c : bundle.hits) {
        uint32_t a_id_o = cluster.particles.getLocalID(hits[c].p_a->id);
        uint32_t b_id_o = cluster.particles.getLocalID(hits[c].p_b->id);

        if (a_id_o != a_id || b_id_o != b_id) {
            throw std::runtime_error("Mismatch ids");
        }

        const Eigen::Matrix4d a_FState_transformation = FStates[a_id].getTransformation();
        const Eigen::Matrix4d b_FState_transformation = FStates[b_id].getTransformation();
        const Eigen::Vector3d a_FState_translation = FStates[a_id].getTranslation();
        const Eigen::Vector3d b_FState_translation = FStates[b_id].getTranslation();

        const Eigen::Vector3d n = contacts[c].global_normal.norm() > 1e-12 ? contacts[c].global_normal.normalized() : (b_FState_translation - a_FState_translation).normalized();

        const double outer_search = hits[c].eps_a + hits[c].eps_b;
        const double outer = outer_search * 0.25;

        const Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, a_FState_transformation);
        const Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, b_FState_transformation);
        const double d1 = p1_B.dot(n) - p1_A.dot(n);

        if (!(d1 > 0.0)) {
            all_d1_above = false;
        }
        // all_d1_above &= d1 > 0.0;

        // double t0 = cluster.particles[a_id].current_state.getTime();
        // double t1 = FStates[a_id].getTime();
        // double ts = t1 - t0;

        // State FState_extrapolate_A = FStates[a_id].extrapolate(ts, hits[c].p_a->getInverseInertiaMatrix());
        // State FState_extrapolate_B = FStates[b_id].extrapolate(ts, hits[c].p_b->getInverseInertiaMatrix());

        // Eigen::Vector3d pf_A = common::transform(contacts[c].local_A, FState_extrapolate_A.getTransformation());
        // Eigen::Vector3d pf_B = common::transform(contacts[c].local_B, FState_extrapolate_B.getTransformation());

        // Eigen::Vector3d v1_A_ext = FStates[a_id].pointVelocity(p1_A, FState_extrapolate_A);
        // Eigen::Vector3d v1_B_ext = FStates[b_id].pointVelocity(p1_B, FState_extrapolate_B);

        const Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, a_inverse_inertia_matrix);
        const Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, b_inverse_inertia_matrix);
        
        const double v1 = v1_A.dot(n) - v1_B.dot(n);

        // Eigen::Vector3d v1_AB_ext = v1_A_ext - v1_B_ext;
        // const double v1 = v1_A_ext.dot(n) - v1_B_ext.dot(n);

        const double m = contacts[c].mass_eff_normal;

        // ======== CANCEL RELATIVE VELOCITY ALONG NORMAL ========

        const double delta_damp = common::lerp(0.5, 0.1, (double)it/(double)max_iterations);
        // v1 = im / m to cancel velocity
        const double delta_impulse_mag = d1 < outer ? v1 * m * delta_damp : -contacts[c].impulse;

        const double last_impulse_mag = contacts[c].impulse;
        contacts[c].impulse = std::max(contacts[c].impulse + delta_impulse_mag, 0.0);

        Eigen::Vector3d total_impulse = n * contacts[c].impulse;
        // globals::logger.printf(3, "Contact %i Total impulse: %f, %f, %f\n", c, total_impulse.x(), total_impulse.y(), total_impulse.z());

        if (it > 0) {
            Eigen::Vector3d average = bundle.normal_average / bundle.hits.size();
            total_impulse = common::lerp(average, total_impulse, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
            total_impulse = d1 < outer ? n.dot(total_impulse) * n : Eigen::Vector3d({0, 0, 0});  // Keep direction of impulse
        }
        else {
            bundle.normal_average += n * last_impulse_mag;
        }
        
        contacts[c].impulse = std::max(total_impulse.dot(n), 0.0);

        const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;

        Eigen::Vector3d update_impulse = n * update_impulse_mag;
        bundle.normal_average += update_impulse;

        if (std::abs(update_impulse_mag) > 1e-8) {
            converged = false;
        }
        update_impulse *= -1;  //  The impulse is calculated A->B but applied B->A so flip it


        const Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, a_FState_translation);
        const Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, b_FState_translation);

        // ========  OFFSET ALONG NORMAL ========

        // const double lower_slop = std::max(outer * 0.25, outer - contacts[c].start_dist);
        // const double slop = common::lerp(lower_slop, outer * 0.25, 0.3);
        const double slop = outer * 0.5;
        const double penetration_depth = outer - d1;

        // d1 = i/m * t offset distance with impulse.  m * d1 / t = i
        const double delta_impulse_offset_mag = d1 < outer ? delta_damp * m * std::max(0.0, penetration_depth - slop) / std::max(t, 1e-6) : -contacts[c].impulse_offset;

        const double last_impulse_offset_mag = contacts[c].impulse_offset;
        contacts[c].impulse_offset = std::max(contacts[c].impulse_offset + delta_impulse_offset_mag, 0.0);
        const double update_impulse_offset_mag = contacts[c].impulse_offset - last_impulse_offset_mag;

        Eigen::Vector3d update_impulse_offset = -n * update_impulse_offset_mag;

        Eigen::Vector3d update_rotational_impulse_offset_A = model::calcTorque(update_impulse_offset, p1_A, a_FState_translation);
        Eigen::Vector3d update_rotational_impulse_offset_B = model::calcTorque(Eigen::Vector3d(-update_impulse_offset), p1_B, b_FState_translation);

        // ======== INTEGRATE ========
        if (!cluster.particles[a_id].is_static)
        {
            impulses[a_id] += update_impulse;
            rotational_impulses[a_id] += update_rotational_impulse_A;
            impulse_offsets[a_id] += update_impulse_offset;
            rotational_impulse_offsets[a_id] += update_rotational_impulse_offset_A;

            const Eigen::Vector3d impulse_a = impulses[a_id];
            const Eigen::Vector3d rotational_impulse_a = rotational_impulses[a_id];
            const Eigen::Vector3d impulse_offset_a = impulse_offsets[a_id];
            const Eigen::Vector3d rotational_impulse_offset_a = rotational_impulse_offsets[a_id];

            FStates[a_id] = updateState(cluster.particles[a_id].future_state, t, impulse_a, impulse_offset_a, rotational_impulse_a, rotational_impulse_offset_a, &cluster.particles[a_id]);
        }
        if (!cluster.particles[b_id].is_static)
        {
            impulses[b_id] -= update_impulse;
            rotational_impulses[b_id] += update_rotational_impulse_B;
            impulse_offsets[b_id] -= update_impulse_offset;
            rotational_impulse_offsets[b_id] += update_rotational_impulse_offset_B;

            const Eigen::Vector3d impulse_b = impulses[b_id];
            const Eigen::Vector3d rotational_impulse_b = rotational_impulses[b_id];
            const Eigen::Vector3d impulse_offset_b = impulse_offsets[b_id];
            const Eigen::Vector3d rotational_impulse_offset_b = rotational_impulse_offsets[b_id];

            FStates[b_id] = updateState(cluster.particles[b_id].future_state, t, impulse_b, impulse_offset_b, rotational_impulse_b, rotational_impulse_offset_b, &cluster.particles[b_id]);
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
                    std::vector<State>& FStates,
                    std::vector<std::vector<collision::ContactBundle>>& colours,
                    std::vector<collision::ContactBundle>& bundles) {
    __itt_string_handle* solve_contacts_task = __itt_string_handle_create("Solve contacts");
    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, solve_contacts_task);

    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    for (collision::ContactBundle& bundle : bundles) {
        bundle.normal_average = {0, 0, 0};
    }
    for (std::vector<collision::ContactBundle>& colour_bundles : colours) {
        for (collision::ContactBundle& bundle : colour_bundles) {
            bundle.normal_average = {0, 0, 0};
        }
    }

    const int max_iterations = globals::opt.sequential_impulse_total_iterations;

    const int parallel_threshold = globals::opt.sequential_parallel_threshold;
    const int parallel_individual_colour_threshold = globals::opt.sequential_parallel_individual_colour_threshold;
    const int parallel_grain_size = globals::opt.sequential_parallel_grain_size;

    std::atomic_bool all_d1_above = true;
    std::atomic_bool converged = false;
    int it = 0;
    while (!converged && it < max_iterations) {
        converged = true;
        all_d1_above = true;
        
        if (colours.size() > 0) {
            for (int colour = 0; colour < colours.size(); colour++) {
                // globals::logger.printf(1, "Colour size: %i\n", colours[colour].size());
                if (colours[colour].size() > parallel_individual_colour_threshold) {
                    // time_point<Clock> start = Clock::now();
                    tbb::parallel_for(tbb::blocked_range<int>(0, colours[colour].size(), parallel_grain_size),
                                    [&](tbb::blocked_range<int> r) {
                        for (int i = r.begin(); i < r.end(); i++) {
                            contact_iteration(colours[colour][i], it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, all_d1_above, converged);
                        }
                    });
                    // time_point<Clock> end = Clock::now();
                    // microseconds diff = duration_cast<microseconds>(end - start);
                    // globals::logger.printf(1, "Parallel size: %i in %iμs\n", colours[colour].size(), diff.count());
                }
                else {
                    // time_point<Clock> start = Clock::now();
                    for (int b = 0; b < colours[colour].size(); b++) {
                        contact_iteration(colours[colour][b], it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, all_d1_above, converged);
                    }
                    // time_point<Clock> end = Clock::now();
                    // microseconds diff = duration_cast<microseconds>(end - start);
                    // globals::logger.printf(1, "Coloured but too small: %i in %iμs\n", colours[colour].size(), diff.count());
                }
            }
        }
        else {
            // time_point<Clock> start = Clock::now();
            for (collision::ContactBundle& bundle : bundles) {
                contact_iteration(bundle, it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, all_d1_above, converged);
            }
            // time_point<Clock> end = Clock::now();
            // microseconds diff = duration_cast<microseconds>(end - start);
            // if (hits.size() > parallel_individual_colour_threshold) {
            //     globals::logger.printf(1, "Not coloured: %i in %iμs\n", hits.size(), diff.count());
            // }
        }

        if (!all_d1_above && converged) {
            globals::logger.printf(4, "%i, Not all d1 above %i\n", cluster_id, it);
        }
        if (!all_d1_above) {
            converged = false;
        }

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

    __itt_task_end(globals::itt_handles.detailed_domain);
    return true;
}

void friction_iteration(collision::ContactBundle& bundle,
                        int it,
                        collision::Cluster& cluster,
                        std::vector<collision::Contact<double>>& hits,
                        std::vector<LocalContact>& contacts,
                        std::vector<Eigen::Vector3d>& impulses,
                        std::vector<Eigen::Vector3d>& rotational_impulses,
                        std::vector<Eigen::Vector3d>& impulse_offsets,
                        std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                        std::vector<State>& FStates,
                        std::atomic_bool& converged) {
    const int max_iterations = globals::opt.sequential_impulse_total_iterations;
    const double t = cluster.step_size;

    std::vector<Eigen::Vector3d> updated_friction_total;
    updated_friction_total.resize(bundle.hits.size());
    for (int c = 0; c < bundle.hits.size(); c++) {
        updated_friction_total[c] = {0, 0, 0};
    }

    const uint32_t a_id = bundle.lower;
    const uint32_t b_id = bundle.upper;

    int i = -1;
    for (int c : bundle.hits) {
        i++;

        const Eigen::Vector3d a_FState_translation = FStates[a_id].getTranslation();
        const Eigen::Vector3d b_FState_translation = FStates[b_id].getTranslation();
        
        const Eigen::Vector3d n = contacts[c].global_normal.norm() > 1e-12 ? contacts[c].global_normal.normalized() : (b_FState_translation - a_FState_translation).normalized();

        const Eigen::Vector3d p0_A = common::transform(contacts[c].local_A, cluster.particles[a_id].current_state.getTransformation());
        const Eigen::Vector3d p0_B = common::transform(contacts[c].local_B, cluster.particles[b_id].current_state.getTransformation());

        const Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
        const Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

        const Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
        const Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());
        const Eigen::Vector3d v1_AB = v1_A - v1_B;
        const Eigen::Vector3d v_tangent_pt = v1_AB - v1_AB.dot(n) * n;
        
        // double t0 = cluster.particles[a_id].current_state.getTime();
        // double t1 = FStates[a_id].getTime();
        // double ts = t1 - t0;

        // State FState_extrapolate_A = FStates[a_id].extrapolate(ts, hits[c].p_a->getInverseInertiaMatrix());
        // State FState_extrapolate_B = FStates[b_id].extrapolate(ts, hits[c].p_b->getInverseInertiaMatrix());

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

        const Eigen::Vector3d v_tangent = v_tangent_pt;
        const Eigen::Vector3d tangent = v_tangent.normalized();

        // if (it == 0) {
        //     globals::logger.printf(3, "%i: %i v_tangent: %f, %f, %f (normal %f, tangent %f)\n", cluster_id, c, v_tangent.x(), v_tangent.y(), v_tangent.z(), contacts[c].impulse, contacts[c].last_friction_impulse.norm());
        //     globals::logger.printf(3, "%i: %i v1_AB: %f, %f, %f first friction\n", cluster_id, c, v1_AB.x(), v1_AB.y(), v1_AB.z());
        //     globals::logger.printf(3, "%i: %i n: %f, %f, %f first friction\n", cluster_id, c, n.x(), n.y(), n.z());
        // }
        
        const double mu = std::min(cluster.particles[a_id].friction_coeff, cluster.particles[b_id].friction_coeff);
        const double max_impulse = (contacts[c].impulse + contacts[c].impulse_offset) * mu;

        double k_tangent = 0.0;
        if (!hits[c].p_a->is_static) {
            const Eigen::Matrix3d ii_A = hits[c].p_a->getInverseInertiaMatrix();
            const Eigen::Vector3d r_A = hits[c].A - hits[c].p_a->future_state.getTranslation();
            const Eigen::Vector3d rt_A = r_A.cross(tangent);
            const double im_A = 1.0 / hits[c].p_a->getMass();
            k_tangent += im_A + rt_A.transpose() * ii_A * rt_A;

        }
        if (!hits[c].p_b->is_static) {
            const Eigen::Matrix3d ii_B = hits[c].p_b->getInverseInertiaMatrix();
            const Eigen::Vector3d r_B = hits[c].B - hits[c].p_b->future_state.getTranslation();
            const Eigen::Vector3d rt_B = r_B.cross(tangent);
            const double im_B = 1.0 / hits[c].p_b->getMass();
            k_tangent += im_B + rt_B.transpose() * ii_B * rt_B;
        }

        // if (k_tangent < 1e-6 && max_impulse > 0) {
        //     globals::logger.printf(1, "Extremely small k_tangent for non-zero impulse");
        // }
        const double m_inv_effective_tangent = k_tangent < 1e-6 ? 0.0 : 1.0 / k_tangent;

        const double friction_damp = common::lerp(0.1, 0.01, (double)it/(double)max_iterations);
        const Eigen::Vector3d delta_friction_impulse = -v_tangent * m_inv_effective_tangent * friction_damp;

        Eigen::Vector3d friction_impulse_total = contacts[c].last_friction_impulse + delta_friction_impulse;
        const Eigen::Vector3d friction_impulse_total_pre_avg = friction_impulse_total;

        if (it > 0) {
            const auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));
            const Eigen::Vector3d average = bundle.tangent_average / bundle.hits.size();
            friction_impulse_total = common::lerp(average, friction_impulse_total, std::min(1.0, std::min(1.0, 2.0*(double)it/(double)max_iterations)));
            friction_impulse_total = -tangent * std::max(0.0, (-tangent).dot(friction_impulse_total));
        }
        friction_impulse_total = friction_impulse_total_pre_avg.normalized() * std::clamp(friction_impulse_total.norm(), 0.0, max_impulse);

        updated_friction_total[i] = friction_impulse_total;

        if ((contacts[c].last_friction_impulse - updated_friction_total[i]).norm() > 1e-8) {
            converged = false;
        }
    }

    bundle.tangent_average = {0, 0, 0};

    i = -1;
    for (int c : bundle.hits) {
        i++;

        bundle.tangent_average += updated_friction_total[i];

        const Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
        const Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

        const Eigen::Vector3d update_friction_impulse = (updated_friction_total[i] - contacts[c].last_friction_impulse) / bundle.hits.size();
        contacts[c].last_friction_impulse = updated_friction_total[i];
        impulses[a_id] += update_friction_impulse;
        impulses[b_id] -= update_friction_impulse;

        const Eigen::Vector3d rotational_impulse_friction_A = model::calcTorque(contacts[c].last_friction_impulse, p1_A, FStates[a_id].getTranslation());
        const Eigen::Vector3d rotational_impulse_friction_B = model::calcTorque(Eigen::Vector3d(-contacts[c].last_friction_impulse), p1_B, FStates[b_id].getTranslation());

        rotational_impulses[a_id] -= contacts[c].last_friction_rotational_impulse_A;
        rotational_impulses[b_id] -= contacts[c].last_friction_rotational_impulse_B;
        contacts[c].last_friction_rotational_impulse_A = rotational_impulse_friction_A;
        contacts[c].last_friction_rotational_impulse_B = rotational_impulse_friction_B;
        rotational_impulses[a_id] += contacts[c].last_friction_rotational_impulse_A;
        rotational_impulses[b_id] += contacts[c].last_friction_rotational_impulse_B;
    }

    if (!cluster.particles[a_id].is_static)
    {
        const Eigen::Vector3d impulse_a = impulses[a_id];
        const Eigen::Vector3d rotational_impulse_a = rotational_impulses[a_id];
        const Eigen::Vector3d impulse_offset_a = impulse_offsets[a_id];
        const Eigen::Vector3d rotational_impulse_offset_a = rotational_impulse_offsets[a_id];

        FStates[a_id] = updateState(cluster.particles[a_id].future_state, t, impulse_a, impulse_offset_a, rotational_impulse_a, rotational_impulse_offset_a, &cluster.particles[a_id]);
    }
    if (!cluster.particles[b_id].is_static)
    {
        const Eigen::Vector3d impulse_b = impulses[b_id];
        const Eigen::Vector3d rotational_impulse_b = rotational_impulses[b_id];
        const Eigen::Vector3d impulse_offset_b = impulse_offsets[b_id];
        const Eigen::Vector3d rotational_impulse_offset_b = rotational_impulse_offsets[b_id];

        FStates[b_id] = updateState(cluster.particles[b_id].future_state, t, impulse_b, impulse_offset_b, rotational_impulse_b, rotational_impulse_offset_b, &cluster.particles[b_id]);
    }
}

bool solve_friction(collision::Cluster& cluster,
                    std::vector<collision::Contact<double>>& hits,
                    std::vector<LocalContact>& contacts,
                    std::vector<Eigen::Vector3d>& impulses,
                    std::vector<Eigen::Vector3d>& rotational_impulses,
                    std::vector<Eigen::Vector3d>& impulse_offsets,
                    std::vector<Eigen::Vector3d>& rotational_impulse_offsets,
                    std::vector<State>& FStates,
                    std::vector<std::vector<collision::ContactBundle>>& colours,
                    std::vector<collision::ContactBundle>& bundles) {
    __itt_string_handle* solve_contacts_task = __itt_string_handle_create("Solve friction");
    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, solve_contacts_task);

    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }
    
    for (collision::ContactBundle& bundle : bundles) {
        bundle.tangent_average = {0, 0, 0};
    }
    for (std::vector<collision::ContactBundle>& colour_bundles : colours) {
        for (collision::ContactBundle& bundle : colour_bundles) {
            bundle.tangent_average = {0, 0, 0};
        }
    }

    const int max_iterations = globals::opt.sequential_impulse_total_iterations;

    const int parallel_threshold = globals::opt.sequential_parallel_threshold;
    const int parallel_individual_colour_threshold = globals::opt.sequential_parallel_individual_colour_threshold;
    const int parallel_grain_size = globals::opt.sequential_parallel_grain_size;

    std::vector<Eigen::Vector3d> tangents;

    std::map<std::pair<int, int>, std::pair<int, Eigen::Vector3d>> average_contact_per_pair;
    
    std::vector<Eigen::Vector3d> updated_friction_total;
    for (int c = 0; c < hits.size(); c++) {
        updated_friction_total.push_back({0, 0, 0});
    }

    std::atomic_bool converged = false;
    int it = 0;
    while (!converged && it < max_iterations || it < 10) {
        converged = true;

        // for (collision::ContactBundle& bundle : bundles) {
        //     friction_iteration(bundle, it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, converged);
        // }

        if (colours.size() > 0) {
            for (int colour = 0; colour < colours.size(); colour++) {
                // globals::logger.printf(1, "Colour size: %i\n", colours[colour].size());
                if (colours[colour].size() > parallel_individual_colour_threshold) {
                    // time_point<Clock> start = Clock::now();
                    tbb::parallel_for(tbb::blocked_range<int>(0, colours[colour].size(), parallel_grain_size),
                                    [&](tbb::blocked_range<int> r) {
                        for (int i = r.begin(); i < r.end(); i++) {
                            friction_iteration(colours[colour][i], it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, converged);
                            // friction_iteration(colours[colour][i], it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, converged);
                        }
                    });
                    // time_point<Clock> end = Clock::now();
                    // microseconds diff = duration_cast<microseconds>(end - start);
                    // globals::logger.printf(1, "Parallel size: %i in %iμs\n", colours[colour].size(), diff.count());
                }
                else {
                    // time_point<Clock> start = Clock::now();
                    for (int b = 0; b < colours[colour].size(); b++) {
                        friction_iteration(colours[colour][b], it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, converged);
                    }
                    // time_point<Clock> end = Clock::now();
                    // microseconds diff = duration_cast<microseconds>(end - start);
                    // globals::logger.printf(1, "Coloured but too small: %i in %iμs\n", colours[colour].size(), diff.count());
                }
            }
        }
        else {
            // time_point<Clock> start = Clock::now();
            for (collision::ContactBundle& bundle : bundles) {
                friction_iteration(bundle, it, cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, converged);
            }
            // time_point<Clock> end = Clock::now();
            // microseconds diff = duration_cast<microseconds>(end - start);
            // if (hits.size() > parallel_individual_colour_threshold) {
            //     globals::logger.printf(1, "Not coloured: %i in %iμs\n", hits.size(), diff.count());
            // }
        }

        it++;
    }

    // if (it == max_iterations) {
    //      __itt_task_end(globals::itt_handles.detailed_domain);
    //     return false;  // TODO this used to cause some very odd errors in the tiny_tower scenario. Still the case?
    // }
    
    globals::logger.printf(3, "%i: Iterations used friction phase: %i\n", cluster_id, it);

    for (int c = 0; c < hits.size(); c++) {
        globals::logger.printf(3, "Contact %i last_impulse_friction: %f, %f, %f (%f)\n", c, contacts[c].last_friction_impulse.x(), contacts[c].last_friction_impulse.y(), contacts[c].last_friction_impulse.z(), contacts[c].last_friction_impulse.norm());
    }
    __itt_task_end(globals::itt_handles.detailed_domain);
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
        std::vector<std::vector<collision::ContactBundle>> colours;
        std::vector<collision::ContactBundle> bundles;
        if (hits.size() > globals::opt.sequential_parallel_threshold) {
            colours = colour_hits(cluster, hits);
        }
        else {
            bundles = bundle_hits(cluster, hits);
        }

        if (!solve_contacts(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, colours, bundles)) return false;

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
            // if (!solve_friction(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, colours, bundles)) return false;
            if (!solve_friction(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, colours, bundles)) return false;

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
            
            if (!solve_contacts(cluster, hits, contacts, impulses, rotational_impulses, impulse_offsets, rotational_impulse_offsets, FStates, colours, bundles)) return false;
            
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

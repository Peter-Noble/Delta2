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

bool SequentialImpulses::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits, bool allow_fail) {
    // __itt_string_handle* sequential_impulse_solve_task = __itt_string_handle_create("Sequential Impulse solve");
    // __itt_string_handle* sequential_impulse_impulses_task = __itt_string_handle_create("Sequential Impulse impulses");
    // __itt_string_handle* sequential_impulse_friction_task = __itt_string_handle_create("Sequential Impulse friction");
    // __itt_string_handle* sequential_impulse_iteration_task = __itt_string_handle_create("Sequential Impulse iteration");
    // __itt_string_handle* sequential_impulse_last_iteration_task = __itt_string_handle_create("Sequential Impulse last iteration");

    // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_solve_task);

    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    std::vector<LocalContact> contacts;

    std::sort(hits.begin(), hits.end());

    for (int c = 0; c < hits.size(); c++) {
        LocalContact lc;

        lc.start_dist = (hits[c].A - hits[c].B).norm();

        if (!(lc.start_dist > 0.0) && cluster.step_size > 1e-4 && allow_fail) {
            // __itt_task_end(globals::itt_handles.detailed_domain);
            if (cluster.step_size < 1e-6) {
                globals::logger.printf(3, "Failed near zero impulse solve\n");
            }
            return false;
        }

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
        lc.last_friction_rotational_impulse_B = { 0.0, 0.0, 0.0 };

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
        if (FStates[p_i].getTransformation().hasNaN()) {
            throw std::runtime_error("Beginning of FStates has NaN");
        }
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
        const double outer_search = hits[c].eps_a + hits[c].eps_b;
        const double outer = outer_search * 0.25;

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

    // Make a first prediction at the future state
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

    std::shared_mutex particle_access_lock;

    // #define USETBB
    #ifndef USETBB
    const int max_iterations = globals::opt.sequential_impulse_total_iterations;
    // globals::logger.printf(3, "Max iterations: %i\n", max_iterations);
    #else
    const int max_iterations = globals::opt.sequential_impulse_total_iterations / globals::opt.sequential_impulse_inner_iterations;
    const int inner_iterations = globals::opt.sequential_impulse_inner_iterations;
    const int grain_size = globals::opt.sequential_impulse_grain_size;
    #endif

    // globals::logger.printf(3, "Max iters: %i\n", max_iterations);

    // globals::logger.printf(3, "hits size: %i, max: %i, inner: %i, grain: %i\n", hits.size(), max_iterations, inner_iterations, grain_size);

    #ifdef USETBB
    Spinlock* mxs[cluster.particles.size()];
    for (int p = 0; p < cluster.particles.size(); p++) {
        mxs[p] = new Spinlock;
    }
    #endif

    if (hits.size() > 0) {
        bool converged = false;
        int it = 0;
        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_impulses_task);
        while (!converged && it < max_iterations) {
            bool all_d1_above = true;
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
            
            #ifndef USETBB
            {
                    for (int c = 0; c < hits.size(); c++) {
            #else
            // tbb::parallel_for(tbb::blocked_range<int>(0, hits.size(), grain_size),
            //     [&](const tbb::blocked_range<int>& r) {
            //         for (int inner_i = 0; inner_i < inner_iterations; inner_i++) {
            //         for (int c = r.begin(); c < r.end(); c++) {
            tbb::task_group task_group;
            task_group.run([&] {
                for (int inner_i = 0; inner_i < inner_iterations; inner_i++) {
                    for (int c = 0; c < hits.size(); c++) {

            #endif
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
                        const double outer = outer_search * 0.25 ;

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
                        const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;

                        Eigen::Vector3d update_impulse = -n * update_impulse_mag;


                        Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, a_FState_translation);
                        Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, b_FState_translation);

                        // ========  OFFSET ALONG NORMAL ========

                        const double lower_slop = std::max(outer * 0.25, outer - contacts[c].start_dist);
                        const double slop = common::lerp(lower_slop, outer * 0.25, 0.3);
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
                }
            #ifdef USETBB
                }
            );
            task_group.wait();
            #endif

            converged = all_d1_above;
            if (converged) {
                for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                    Eigen::Vector3d impulse_diff = impulses[p_i] - last_impulse[p_i];
                    if (impulse_diff.norm() > 1e-8) {
                        converged = false;
                        break;
                    }
                }
            }

            for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                Eigen::Vector3d impulse_diff = impulses[p_i] - last_impulse[p_i];
                // globals::logger.printf(3, "impulse: %f, %f, %f\n", impulses[p_i].x(), impulses[p_i].y(), impulses[p_i].z());
                // globals::logger.printf(3, "impulse diff: %f, %f, %f\n", impulse_diff.x(), impulse_diff.y(), impulse_diff.z());
            }

            it++;
            // __itt_task_end(globals::itt_handles.detailed_domain);
        }

        globals::logger.printf(3, "%i: Iterations used: %i\n", cluster_id, it);

        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            assert(!impulses[p_i].hasNaN());
            assert(impulses[p_i].norm() < 10000);
            if (impulses[p_i].norm() > 10000) {
                globals::logger.printf(3, "Really big impulse\n");
            }
            if (impulses[p_i].hasNaN()) {
                // // throw std::runtime_error("Nan!");
                globals::logger.printf(3, "Nan!\n");
                return false;
            }
        }

        // __itt_task_end(globals::itt_handles.detailed_domain);

        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_friction_task);

        std::vector<Eigen::Vector3d> updated_friction_total;
        for (int c = 0; c < hits.size(); c++) {
            updated_friction_total.push_back({0, 0, 0});
        }

        for (int it = 0; it < 400; it++) {
            std::map<std::pair<int, int>, int> contacts_per_pair;
            // globals::logger.printf(3, "===============================\n");

            std::vector<common::Edge<double>> view_hits;
            bool show = false;

            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                Eigen::Vector3d n = contacts[c].global_normal.normalized();
                if (n.hasNaN()) {
                    Eigen::Vector3d a_FState_translation = FStates[a_id].getTranslation();
                    Eigen::Vector3d b_FState_translation = FStates[b_id].getTranslation();
                    n = (b_FState_translation - a_FState_translation).normalized();
                }

                const double outer_search = hits[c].eps_a + hits[c].eps_b;
                const double outer = outer_search * 0.25;

                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());
                const double d1 = p1_B.dot(n) - p1_A.dot(n);

                Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, hits[c].p_a->getInverseInertiaMatrix());
                Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, hits[c].p_b->getInverseInertiaMatrix());

                const double m = contacts[c].mass_eff_normal;

                Eigen::Vector3d v1_AB = v1_A - v1_B;

                Eigen::Vector3d v_tangent = v1_AB - v1_AB.dot(n) * n;
                double mu = std::min(cluster.particles[a_id].friction_coeff, cluster.particles[b_id].friction_coeff);
                double max_force = contacts[c].impulse * mu * 1000;

                // globals::logger.printf(3, "it: %i, c: %i, max_force: %.4f, v_tangent: %.4f\n", it, c, max_force, v_tangent.norm());
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

                Eigen::Vector3d delta_friction_impulse = -v_tangent * m_inv_effective_tangent * 0.25;
                // if (max_force > 1e-6) {
                //     globals::logger.printf(3, "delta_friction_impulse: %.4f, %.4f, %.4f\n", delta_friction_impulse.x(), delta_friction_impulse.y(), delta_friction_impulse.z());
                // }
                Eigen::Vector3d friction_impulse_total = contacts[c].last_friction_impulse + delta_friction_impulse;
                friction_impulse_total = friction_impulse_total.normalized() * std::clamp(friction_impulse_total.norm(), 0.0, max_force);

                updated_friction_total[c] = friction_impulse_total;

                if (updated_friction_total[c].hasNaN()) {
                    globals::logger.printf(3, "c: %i\n", c);
                    globals::logger.printf(3, FStates[a_id].getTransformation().hasNaN() ? "FStates a_id has nan" : "FStates a_id is valid\n");
                    globals::logger.printf(3, FStates[b_id].getTransformation().hasNaN() ? "FStates b_id has nan" : "FStates b_id is valid\n");
                    globals::logger.printf(3, "contacts[c].local_A: %f, %f, %f\n", contacts[c].local_A.x(), contacts[c].local_A.y(), contacts[c].local_A.z());
                    globals::logger.printf(3, "contacts[c].local_B: %f, %f, %f\n", contacts[c].local_B.x(), contacts[c].local_B.y(), contacts[c].local_B.z());
                    globals::logger.printf(3, "p1_A: %f, %f, %f\n", p1_A.x(), p1_A.y(), p1_A.z());
                    globals::logger.printf(3, "p1_B: %f, %f, %f\n", p1_B.x(), p1_B.y(), p1_B.z());
                    globals::logger.printf(3, "v1_AB: %f, %f, %f\n", v1_AB.x(), v1_AB.y(), v1_AB.z());
                    globals::logger.printf(3, "rt_A: %f, %f, %f\n", rt_A.x(), rt_A.y(), rt_A.z());
                    globals::logger.printf(3, "rt_B: %f, %f, %f\n", rt_B.x(), rt_B.y(), rt_B.z());
                    globals::logger.printf(3, "tangent: %f, %f, %f\n", tangent.x(), tangent.y(), tangent.z());
                    globals::logger.printf(3, "k_tangent: %f\n", k_tangent);
                    globals::logger.printf(3, "m_inv_effective_tangent: %f\n", m_inv_effective_tangent);
                    globals::logger.printf(3, "delta_friction_impulse: %f, %f, %f\n", delta_friction_impulse.x(), delta_friction_impulse.y(), delta_friction_impulse.z());
                    // throw std::runtime_error("Nans");
                }

                auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));
                if (contacts_per_pair.find(key) != contacts_per_pair.end()) {
                    contacts_per_pair[key]++;
                }
                else {
                    contacts_per_pair[key] = 1;
                }

                // if (it == 99 && max_force > 1e-6) {
                //     globals::logger.printf(3, "v_tangent: %.4f, fric impulse: %.4f, max_force: %.4f\n", v_tangent.norm(), friction_impulse_total.norm(), max_force);
                //     show |= v_tangent.norm() > 1e-3;
                //     view_hits.push_back(common::Edge(hits[c].A, hits[c].B));                    
                // }
            }

            // if (show) {
            //     common::Viewer view;
            //     view.addParticleFuture(cluster.particles[0]);
            //     view.addParticleFuture(cluster.particles[1]);
            //     for (auto e : view_hits) {
            //         view.addEdge(e);
            //     }    
            //     view.show();
            // }

            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
                auto key = std::make_pair(std::min(a_id, b_id), std::max(a_id, b_id));

                Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, FStates[a_id].getTransformation());
                Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, FStates[b_id].getTransformation());

                Eigen::Vector3d update_friction_impulse = (updated_friction_total[c] - contacts[c].last_friction_impulse) / contacts_per_pair[key];
                contacts[c].last_friction_impulse = updated_friction_total[c];
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
            }

            // ======== INTEGRATE ========
            for (int p = 0; p < cluster.particles.size(); p++) {
                if (!cluster.particles[p].is_static)
                {
                    Eigen::Vector3d impulse = impulses[p];
                    Eigen::Vector3d rotational_impulse = rotational_impulses[p];
                    Eigen::Vector3d impulse_offset = impulse_offsets[p];
                    Eigen::Vector3d rotational_impulse_offset = rotational_impulse_offsets[p];
                    FStates[p] = updateState(cluster.particles[p].future_state, t, impulse, impulse_offset, rotational_impulse, rotational_impulse_offset, &cluster.particles[p]);
                }
            }
        }

        // __itt_task_end(globals::itt_handles.detailed_domain);

        for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
            assert(!impulses[p_i].hasNaN());
            assert(impulses[p_i].norm() < 10000);
            if (impulses[p_i].norm() > 10000) {
                globals::logger.printf(3, "Really big impulse\n");
            }
            if (impulses[p_i].hasNaN()) {
                // throw std::runtime_error("Nan!");
                globals::logger.printf(3, "Nan!\n");
                return false;
            }
        }



        it = 0;
        converged = false;
        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, sequential_impulse_impulses_task);
        while (!converged && it < max_iterations) {
            bool all_d1_above = true;
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
            
            #ifndef USETBB
            {
                    for (int c = 0; c < hits.size(); c++) {
            #else
            // tbb::parallel_for(tbb::blocked_range<int>(0, hits.size(), grain_size),
            //     [&](const tbb::blocked_range<int>& r) {
            //         for (int inner_i = 0; inner_i < inner_iterations; inner_i++) {
            //         for (int c = r.begin(); c < r.end(); c++) {
            tbb::task_group task_group;
            task_group.run([&] {
                for (int inner_i = 0; inner_i < inner_iterations; inner_i++) {
                    for (int c = 0; c < hits.size(); c++) {

            #endif
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
                        const double outer = outer_search * 0.25 ;

                        Eigen::Vector3d p1_A = common::transform(contacts[c].local_A, a_FState_transformation);
                        Eigen::Vector3d p1_B = common::transform(contacts[c].local_B, b_FState_transformation);
                        const double d1 = p1_B.dot(n) - p1_A.dot(n);

                        all_d1_above &= d1 > 0.0;

                        Eigen::Vector3d v1_A = FStates[a_id].pointVelocity(p1_A, a_inverse_inertia_matrix);
                        Eigen::Vector3d v1_B = FStates[b_id].pointVelocity(p1_B, b_inverse_inertia_matrix);
                        const double v1 = v1_A.dot(n) - v1_B.dot(n);

                        const double m = contacts[c].mass_eff_normal;

                        // ======== CANCEL RELATIVE VELOCITY ALONG NORMAL ========

                        const double delta_damp = common::lerp(0.5, 0.01, (double)it/(double)max_iterations);
                        // v1 = im / m to cancel velocity
                        const double delta_impulse_mag = d1 < outer ? v1 * m * delta_damp : -contacts[c].impulse;

                        double last_impulse_mag = contacts[c].impulse;
                        contacts[c].impulse = std::max(contacts[c].impulse + delta_impulse_mag, 0.0);
                        const double update_impulse_mag = contacts[c].impulse - last_impulse_mag;

                        Eigen::Vector3d update_impulse = -n * update_impulse_mag;

                        if (update_impulse.hasNaN()) {
                            globals::logger.printf(2, "Update impulse has nan\n");
                            return false;
                            throw std::runtime_error("Update impulse has nan");
                        }

                        Eigen::Vector3d update_rotational_impulse_A = model::calcTorque(update_impulse, p1_A, a_FState_translation);
                        Eigen::Vector3d update_rotational_impulse_B = model::calcTorque(Eigen::Vector3d(-update_impulse), p1_B, b_FState_translation);

                        // ========  OFFSET ALONG NORMAL ========

                        const double lower_slop = std::max(outer * 0.25, outer - contacts[c].start_dist);
                        const double slop = common::lerp(lower_slop, outer * 0.25, 0.3);
                        const double penetration_depth = outer - d1;

                        // d1 = i/m * t offset distance with impulse.  m * d1 / t = i
                        double delta_impulse_offset_mag = d1 < outer ? delta_damp * m * std::max(0.0, penetration_depth - slop) / std::max(t, 1e-6) : -contacts[c].impulse_offset;
                        // if (t < 1e-8) {
                        //     delta_impulse_offset_mag = 0;
                        // }

                        if (std::isnan(delta_impulse_offset_mag)) {
                            throw std::runtime_error("offset mag is nan");
                        }
                        if (std::isinf(delta_impulse_offset_mag)) {
                            throw std::runtime_error("offset mag is inf");
                        }

                        double last_impulse_offset_mag = contacts[c].impulse_offset;
                        contacts[c].impulse_offset = std::max(contacts[c].impulse_offset + delta_impulse_offset_mag, 0.0);
                        const double update_impulse_offset_mag = contacts[c].impulse_offset - last_impulse_offset_mag;

                        if (std::isinf(update_impulse_offset_mag)) {
                            throw std::runtime_error("update impulse offset mag has inf");
                        }

                        Eigen::Vector3d update_impulse_offset = -n * update_impulse_offset_mag;

                        if (update_impulse_offset.hasNaN()) {
                            throw std::runtime_error("Update impulse offset has nan");
                        }
                        if (a_FState_translation.hasNaN()) {
                            throw std::runtime_error("a_FState_translation has nan");
                        }

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

                                if (impulse_offsets[a_id].hasNaN()) {
                                    throw std::runtime_error("impulse offsets has nan");
                                }

                                if (update_impulse_offset.hasNaN()) {
                                    throw std::runtime_error("Update impulse offset has nan");
                                }

                                if (update_rotational_impulse_offset_A.hasNaN()) {
                                    globals::logger.printf(3, "uio: (%f, %f, %f)\n", update_impulse_offset.x(), update_impulse_offset.y(), update_impulse_offset.z());
                                    globals::logger.printf(3, "p1_A: (%f, %f, %f)\n", p1_A.x(), p1_A.y(), p1_A.z());
                                    globals::logger.printf(3, "aFt: (%f, %f, %f)\n", a_FState_translation.x(), a_FState_translation.y(), a_FState_translation.z());

                                    throw std::runtime_error("Update impulse offset has nan");
                                }

                                impulses[a_id] += update_impulse;
                                rotational_impulses[a_id] += update_rotational_impulse_A;
                                impulse_offsets[a_id] += update_impulse_offset;
                                rotational_impulse_offsets[a_id] += update_rotational_impulse_offset_A;

                                if (impulse_offsets[a_id].hasNaN()) {
                                    throw std::runtime_error("Impulse offset has nan");
                                }
                                if (rotational_impulse_offsets[a_id].hasNaN()) {
                                    globals::logger.printf(3, "uio: (%f, %f, %f)\n", update_impulse_offset.x(), update_impulse_offset.y(), update_impulse_offset.z());
                                    globals::logger.printf(3, "p1_A: (%f, %f, %f)\n", p1_A.x(), p1_A.y(), p1_A.z());
                                    globals::logger.printf(3, "aFt: (%f, %f, %f)\n", a_FState_translation.x(), a_FState_translation.y(), a_FState_translation.z());

                                    throw std::runtime_error("Impulse rotational offset has nan");
                                }

                                Eigen::Vector3d impulse_a = impulses[a_id];
                                Eigen::Vector3d rotational_impulse_a = rotational_impulses[a_id];
                                Eigen::Vector3d impulse_offset_a = impulse_offsets[a_id];
                                Eigen::Vector3d rotational_impulse_offset_a = rotational_impulse_offsets[a_id];

                                if (impulse_a.hasNaN() || rotational_impulse_a.hasNaN()) {
                                    // throw std::runtime_error("Nan in a state update\n");
                                }
                                if (impulse_offset_a.hasNaN()) {
                                    // throw std::runtime_error("Nan in a state update impulse offset\n");
                                }
                                if (rotational_impulse_offset_a.hasNaN()) {
                                    globals::logger.printf(3, "uio: (%f, %f, %f)\n", update_impulse_offset.x(), update_impulse_offset.y(), update_impulse_offset.z());
                                    globals::logger.printf(3, "p1_A: (%f, %f, %f)\n", p1_A.x(), p1_A.y(), p1_A.z());
                                    globals::logger.printf(3, "aFt: (%f, %f, %f)\n", a_FState_translation.x(), a_FState_translation.y(), a_FState_translation.z());

                                    // throw std::runtime_error("Nan in a state update rotational impulse offset\n");
                                }
                                FStates[a_id] = updateState(cluster.particles[a_id].future_state, t, impulse_a, impulse_offset_a, rotational_impulse_a, rotational_impulse_offset_a, &cluster.particles[a_id]);
                                
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

                                if (impulse_b.hasNaN() || rotational_impulse_b.hasNaN() || impulse_offset_b.hasNaN() || rotational_impulse_offset_b.hasNaN()) {
                                    // throw std::runtime_error("Nan in b state update\n");
                                }
                                FStates[b_id] = updateState(cluster.particles[b_id].future_state, t, impulse_b, impulse_offset_b, rotational_impulse_b, rotational_impulse_offset_b, &cluster.particles[b_id]);
                                
                                // #ifdef USETBB
                                // mxs[b_id]->unlock();
                                // #endif
                            }
                        }
                    }
                }
            #ifdef USETBB
                }
            );
            task_group.wait();
            #endif

            converged = all_d1_above;
            if (converged) {
                for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                    Eigen::Vector3d impulse_diff = impulses[p_i] - last_impulse[p_i];
                    if (impulse_diff.norm() > 1e-8) {
                        converged = false;
                        if (it > 500) {
                            globals::logger.printf(3, "Diff %f\n", impulse_diff.norm());
                        }
                        break;
                    }
                }
            }
            for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
                Eigen::Vector3d impulse_diff = impulses[p_i] - last_impulse[p_i];
                if (impulses[p_i].hasNaN()) {
                    throw std::runtime_error("Impulse has NaN\n");
                }
                // globals::logger.printf(3, "impulse: %f, %f, %f\n", impulses[p_i].x(), impulses[p_i].y(), impulses[p_i].z());
                // globals::logger.printf(3, "impulse diff: %f, %f, %f\n", impulse_diff.x(), impulse_diff.y(), impulse_diff.z());
            }

            if (converged && !all_d1_above) {
                globals::logger.printf(3, "Incorrectly converged\n");
                throw std::runtime_error("Incorrectly converged");
            }

            it++;
            // __itt_task_end(globals::itt_handles.detailed_domain);
        }

        globals::logger.printf(3, "%i: Iterations used 2nd phase: %i\n", cluster_id, it);





        {
            std::lock_guard draws_lock(globals::contact_draws_lock);
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);
                // globals::logger.printf(3, "Contact with ids: %i(%i), %i(%i), %f, (%f, %f, %f)\n", a_id, hits[c].p_a->id, b_id, hits[c].p_b->id, contacts[c].impulse, contacts[c].global_normal.x(), contacts[c].global_normal.y(), contacts[c].global_normal.z());
                globals::contact_draws.push_back(std::make_tuple(contacts[c].global_centre, contacts[c].global_normal.normalized() * contacts[c].impulse, contacts[c].last_friction_impulse));
            }
        }
        // globals::logger.printf(3, "Finish contact in %i iterations\n", it);
    }

    #ifdef USETBB
    for (int p = 0; p < cluster.particles.size(); p++) {
        delete mxs[p];
    }
    #endif

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        // assert(!impulses[p_i].hasNaN());
        // assert(impulses[p_i].norm() < 10000);
        if (impulses[p_i].norm() > 10000) {
            globals::logger.printf(3, "Really big impulse\n");
            return false;
        }
        if (impulses[p_i].hasNaN()) {
            // // throw std::runtime_error("Nan!");
            globals::logger.printf(3, "Nan!\n");
            return false;
        }
        if (FStates[p_i].getTransformation().hasNaN()) {
            // throw std::runtime_error("Nan in transformation FState");
        }
    }

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        cluster.particles[p_i].future_state = FStates[p_i];
    }

    _warm_start.update(contacts, hits);

    // common::Viewer view;
    // for (int p = 0; p < cluster.particles.size(); p++) {
    //     view.addParticleInterval(cluster.particles[p]);
    // }
    // for (int c = 0; c < hits.size(); c++) {
    //     view.addEdge(common::Edge(hits[c].A, hits[c].B));
    // }
    // view.show();

    // __itt_task_end(globals::itt_handles.detailed_domain);

    return true;
}

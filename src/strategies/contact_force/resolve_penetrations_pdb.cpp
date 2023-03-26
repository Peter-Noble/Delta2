#include "resolve_penetrations_pbd.h"
#include "../../collision_detection/contact.h"
#include "../../collision_detection/full_tree_comparison.h"
#include "../../common/common.h"
#include "../../strategies/contact_force/LocalContact.h"
#include "../../strategies/contact_detection/contact_detection_comparison.h"
#include "../../model/forces.h"
#include "../contact_force/update_state.h"
#include "../../globals.h"

using namespace Delta2;
using namespace Eigen;

void Delta2::collision::resolvePenetrationsPBD(collision::Cluster cluster, bool future_only) {
    int cluster_id = -1;
    for (Particle* p : cluster.particles) {
        if (!p->is_static) {
            cluster_id = p->cluster_id;
            break;
        }
    }

    globals::logger.printf(2, "%i: Resolving penetrations\n", cluster_id);

    // common::Viewer view;
    // for (Particle* p : cluster.particles) {
    //     view.addParticleInterval(*p);
    // }
    // view.show();

    std::vector<collision::Contact<double>> hits;

    strategy::ContactDetectionComparison contact_detection;

    bool new_full_comparison = true;

    double mult = 1.0;

    int i = 0;

    while (new_full_comparison) {
        globals::logger.printf(3, "%i, New full comparison %i\n", cluster_id, i);
        i++;
        new_full_comparison = false;
        // globals::logger.printf(3, "Resolve full comparison\n");
        for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
        {
            collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

            int a_id = cluster.particles.getLocalID(b.A.first); // geo id
            int b_id = cluster.particles.getLocalID(b.B.first);

            Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
            Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

            std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], contact_detection);
            std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

            for (collision::Contact<double> &c : filtered)
            {
                hits.push_back(c);
            }
        }
        
        std::vector<LocalContact> contacts;

        std::sort(hits.begin(), hits.end());

        for (int c = 0; c < hits.size(); c++) {
            LocalContact lc;

            lc.start_dist = (hits[c].A - hits[c].B).norm();

            Matrix4d A_T_i = hits[c].p_a->future_state.getTransformation().inverse();
            lc.local_A = common::transform(hits[c].A, A_T_i);
            Matrix4d B_T_i = hits[c].p_b->future_state.getTransformation().inverse();
            lc.local_B = common::transform(hits[c].B, B_T_i);

            lc.global_centre = (hits[c].A + hits[c].B) / 2.0;
            lc.global_normal = (hits[c].B - hits[c].A) / 2.0;
            lc.global_tangent = {0, 0, 0}; // todo

            Vector3d r_A = hits[c].A - hits[c].p_a->future_state.getTranslation();
            Vector3d r_B = hits[c].B - hits[c].p_b->future_state.getTranslation();

            Vector3d rn_A = r_A.cross(lc.global_normal);
            Vector3d rn_B = r_B.cross(lc.global_normal);

            Vector3d rt_A = r_A.cross(lc.global_tangent);
            Vector3d rt_B = r_B.cross(lc.global_tangent);

            double im_A = 1.0 / hits[c].p_a->getMass();
            double im_B = 1.0 / hits[c].p_b->getMass();
            Matrix3d ii_A = hits[c].p_a->getInverseInertiaMatrix();
            Matrix3d ii_B = hits[c].p_b->getInverseInertiaMatrix();

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

        int inner = 0;
        bool changed = true;
        while (changed && inner < 100) {
            inner++;
            changed = false;
            // globals::logger.printf(3, "Resolve inner loop\n");
            for (int c = 0; c < hits.size(); c++) {
                uint32_t a_id = cluster.particles.getLocalID(hits[c].p_a->id);
                uint32_t b_id = cluster.particles.getLocalID(hits[c].p_b->id);

                const double outer_search = hits[c].eps_a + hits[c].eps_b;
                const double outer = 1e-4;

                Vector3d p1_A = common::transform(contacts[c].local_A, cluster.particles[a_id].future_state.getTransformation());
                Vector3d p1_B = common::transform(contacts[c].local_B, cluster.particles[b_id].future_state.getTransformation());

                Vector3d n = contacts[c].global_normal.normalized();
                if (n.hasNaN() || n.norm() < 1e-8) {
                    // n = (hits[c].p_b->future_state.getTranslation() - hits[c].p_a->future_state.getTranslation()).normalized();

                    Eigen::Matrix3d a_inverse_inertia_matrix = hits[c].p_a->getInverseInertiaMatrix();
                    Eigen::Matrix3d b_inverse_inertia_matrix = hits[c].p_b->getInverseInertiaMatrix();
                    Eigen::Vector3d v1_A = cluster.particles[a_id].future_state.pointVelocity(p1_A, a_inverse_inertia_matrix);
                    Eigen::Vector3d v1_B = cluster.particles[b_id].future_state.pointVelocity(p1_B, b_inverse_inertia_matrix);

                    if ((v1_A - v1_B).norm() > 1e-8) {
                        n = (v1_A - v1_B).normalized();
                    }
                    else {
                        n = (hits[c].p_b->future_state.getTranslation() - hits[c].p_a->future_state.getTranslation()).normalized();
                    }
                }

                // const double d1 = p1_B.dot(n) - p1_A.dot(n);
                const double d1 = (p1_B - p1_A).norm();

                // globals::logger.printf(3, "p1_A: (%f, %f, %f)\n", p1_A.x(), p1_A.y(), p1_A.z());
                // globals::logger.printf(3, "p1_B: (%f, %f, %f)\n", p1_B.x(), p1_B.y(), p1_B.z());

                // globals::logger.printf(3, "n: (%f, %f, %f)\n", n.x(), n.y(), n.z());

                // globals::logger.printf(3, "d1: %f\n", d1);

                if (d1 < outer) {
                    changed = true;
                    new_full_comparison = true;
                    if (i > 100) {
                        globals::logger.printf(2, "Can't fix positions so scattering\n");
                        if (!cluster.particles[a_id].is_static) {
                            cluster.particles[a_id].current_state.setTranslation(Eigen::Vector3d({globals::opt.rand_float(1000000), globals::opt.rand_float(1000000), globals::opt.rand_float(1000000)}));
                            cluster.particles[a_id].future_state.setTranslation(Eigen::Vector3d({globals::opt.rand_float(1000000), globals::opt.rand_float(1000000), globals::opt.rand_float(1000000)}));
                        }
                        if (!cluster.particles[b_id].is_static) {
                            cluster.particles[b_id].current_state.setTranslation(Eigen::Vector3d({globals::opt.rand_float(1000000), globals::opt.rand_float(1000000), globals::opt.rand_float(1000000)}));
                            cluster.particles[b_id].future_state.setTranslation(Eigen::Vector3d({globals::opt.rand_float(1000000), globals::opt.rand_float(1000000), globals::opt.rand_float(1000000)}));
                        }
                        continue;
                    }

                    const double m = contacts[c].mass_eff_normal;

                    Vector3d a_correct_dir = hits[c].p_a->future_state.getTranslation() - contacts[c].global_centre;
                    Vector3d b_correct_dir = hits[c].p_b->future_state.getTranslation() - contacts[c].global_centre;

                    Vector3d correct_dir = (hits[c].p_a->future_state.getTranslation() - hits[c].p_b->future_state.getTranslation()).normalized();

                    const double penetration_depth = outer - d1;

                    // d1 = i/m * t offset distance with impulse.  m * d1 / t = i
                    const double t = 1;
                    const double delta_impulse_offset_mag = mult * 0.5 * m * std::max(0.0, penetration_depth) / t;

                    // globals::logger.printf(3, "offset delta: %f\n", delta_impulse_offset_mag);

                    Vector3d update_impulse_offset = n * delta_impulse_offset_mag;

                    // globals::logger.printf(3, "uio: (%f, %f, %f)\n", update_impulse_offset.x(), update_impulse_offset.y(), update_impulse_offset.z());

                    Vector3d update_rotational_impulse_offset_A = model::calcTorque(Vector3d(-update_impulse_offset), p1_A, hits[c].p_a->future_state.getTranslation());
                    Vector3d update_rotational_impulse_offset_B = model::calcTorque(update_impulse_offset, p1_B, hits[c].p_b->future_state.getTranslation());

                    if (!cluster.particles[a_id].is_static)
                    {
                        Vector3d aT = cluster.particles[a_id].future_state.getTranslation();
                        // globals::logger.printf(3, "pid pre: %i (%f, %f, %f)\n", a_id, aT.x(), aT.y(), aT.z());

                        Vector3d zero = Vector3d::Zero();
                        cluster.particles[a_id].future_state = updateState(cluster.particles[a_id].future_state, t, zero, -update_impulse_offset, zero, update_rotational_impulse_offset_A, &cluster.particles[a_id]);

                        if (!future_only) {
                            cluster.particles[a_id].current_state = updateState(cluster.particles[a_id].current_state, t/2, zero, update_impulse_offset, zero, update_rotational_impulse_offset_A, &cluster.particles[a_id]);
                            cluster.particles[a_id].last_state = updateState(cluster.particles[a_id].last_state, t/4, zero, update_impulse_offset, zero, update_rotational_impulse_offset_A, &cluster.particles[a_id]);
                        }

                        aT = cluster.particles[a_id].future_state.getTranslation();
                        // globals::logger.printf(3, "pid post: %i (%f, %f, %f)\n", a_id, aT.x(), aT.y(), aT.z());
                    }

                    if (!cluster.particles[b_id].is_static)
                    {
                        Vector3d bT = cluster.particles[b_id].future_state.getTranslation();
                        // globals::logger.printf(3, "pid pre: %i (%f, %f, %f)\n", b_id, bT.x(), bT.y(), bT.z());

                        Vector3d zero = Vector3d::Zero();
                        cluster.particles[b_id].future_state = updateState(cluster.particles[b_id].future_state, t, zero, update_impulse_offset, zero, update_rotational_impulse_offset_B, &cluster.particles[b_id]);

                        if (!future_only) {
                            cluster.particles[b_id].current_state = updateState(cluster.particles[b_id].current_state, t/2, zero, -update_impulse_offset, zero, update_rotational_impulse_offset_B, &cluster.particles[b_id]);
                            cluster.particles[b_id].last_state = updateState(cluster.particles[b_id].last_state, t/4, zero, -update_impulse_offset, zero, update_rotational_impulse_offset_B, &cluster.particles[b_id]);
                        }

                        bT = cluster.particles[b_id].future_state.getTranslation();
                        // globals::logger.printf(3, "pid post: %i (%f, %f, %f)\n", b_id, bT.x(), bT.y(), bT.z());
                    }
                }
            }
        }

        // common::Viewer view;
        // for (int p = 0; p < cluster.particles.size(); p++) {
        //     view.addParticleInterval(cluster.particles[p]);
        // }
        // view.show();

        hits.clear();
        mult *= 1.1;
    }
}

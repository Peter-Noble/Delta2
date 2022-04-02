#include "separate_clusters.h"
#include "../collision_detection/full_tree_comparison.h"
#include "igl/connected_components.h"
#include <chrono>
#include <omp.h>
#include <ittnotify.h>

using namespace Delta2;

/*
 *  toc on the resulting collisions is 
 */
std::vector<collision::ContinuousContact<double>> collision::separateContinuousCollisionClusters(collision::BroadPhaseCollisions& broad_phase, std::vector<Delta2::Particle>& particles, double max_time, std::vector<std::vector<Particle*>>& cluster_particles_out, std::vector<collision::BroadPhaseCollisions>& cluster_interactions_out, std::vector<bool>& sleeping) {
    __itt_domain* domain = __itt_domain_create("My Domain");

    int num_particles = particles.size();
    Eigen::SparseMatrix<int> interaction_graph(num_particles, num_particles);

    std::mutex lock;

    for (int b_i = 0; b_i < broad_phase.size(); b_i++) {
        collision::BroadPhaseCollision& b = broad_phase[b_i];

        int a_id = b.first.first; // geo id
        int b_id = b.second.first;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            // std::lock_guard<std::mutex> guard(lock);
            interaction_graph.coeffRef(a_id, b_id) = 1;
            interaction_graph.coeffRef(b_id, a_id) = 1;
        }
    }

    Eigen::MatrixXi components;
    Eigen::MatrixXi sizes;
    igl::connected_components(interaction_graph, components, sizes);

    cluster_particles_out.clear();
    cluster_interactions_out.clear();
    sleeping.clear();

    for (int s = 0; s < sizes.rows(); s++) {
        cluster_particles_out.emplace_back();
        cluster_particles_out[s].reserve(sizes(s, 0));
        cluster_interactions_out.emplace_back();
        sleeping.push_back(true);
    }

    for (int p_i = 0; p_i < num_particles; p_i++) {
        Particle* address = &(particles[p_i]);
        cluster_particles_out[components(p_i, 0)].push_back(address);
        if (!particles[p_i].getSleeping()) {
            sleeping[components(p_i, 0)] = false;
        }
    }

    #pragma omp parallel
    {
        // for (collision::BroadPhaseCollision& b : broad_phase) {
        #pragma omp single
        {
            #pragma omp taskloop
            for (int b_i = 0; b_i < broad_phase.size(); b_i++) {
                __itt_string_handle* continuous_pair_task = __itt_string_handle_create("Continuous pair");
                __itt_task_begin(domain, __itt_null, __itt_null, continuous_pair_task);

                int tid = omp_get_thread_num();

                collision::BroadPhaseCollision& b = broad_phase[b_i];

                int a_id = b.first.first; // geo id
                int b_id = b.second.first;

                Eigen::Vector3d a_centre = particles[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = particles[b_id].current_state.getTranslation();

                double max_time_step_for_pair = std::min(particles[a_id].last_time_step_size, particles[b_id].last_time_step_size); // Use the minimum of the maximum allowed timesteps

                double min_time;

                __itt_string_handle* continuous_trees_continuous_task = __itt_string_handle_create("Compare trees full continuous");
                __itt_task_begin(domain, __itt_null, __itt_null, continuous_trees_continuous_task);
                std::vector<collision::ContinuousContact<double>> Ccs = collision::compareTreesFullContinuous<double, 8, 8>(particles[a_id], particles[b_id], max_time_step_for_pair, min_time);
                __itt_task_end(domain);

                if (Ccs.size() > 0) {
                    double min_scaling_seen = 1.0;
                    double min_toc = 1.0;
                    for (collision::ContinuousContact<double>& c : Ccs) {
                        
                    }
                    printf("Min toc: %f\n", min_toc);
                }
                __itt_task_end(domain);
            }
        }
    }

    for (collision::BroadPhaseCollision& b : broad_phase) {
        int a_id = b.first.first; // geo id
        int b_id = b.second.first;
        int id;
        int static_id = -1;
        int cluster_id;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            assert(components(a_id, 0) == components(b_id, 0));
            id = a_id;
        }
        else if (particles[a_id].is_static && particles[b_id].is_static) {
            continue;
        }
        else if (particles[a_id].is_static) {
            id = b_id;
            static_id = a_id;
            cluster_id = components(b_id, 0);
        }
        else {
            id = a_id;
            static_id = b_id;
            cluster_id = components(a_id, 0);
        }

        if (static_id >= 0) {
            Particle* address = &(particles[static_id]);
            if (std::find(cluster_particles_out[cluster_id].begin(), cluster_particles_out[cluster_id].end(), address) == cluster_particles_out[cluster_id].end()) {
                cluster_particles_out[cluster_id].push_back(address);
            }
        }

        cluster_interactions_out[components(a_id, 0)].push_back(b);
    }
}

// TODO Split into at least 2 parts.
//  1. Broadphase separate connected components
//  2. CompareTreesContinuous
void collision::separateCollisionClustersWithTimeStepSelection(collision::BroadPhaseCollisions& broad_phase, std::vector<Delta2::Particle>& particles, double max_time, std::vector<std::vector<Particle*>>& cluster_particles_out, std::vector<collision::BroadPhaseCollisions>& cluster_interactions_out, std::vector<double>& cluster_step_size_out, std::vector<bool>& sleeping) {
    __itt_domain* domain = __itt_domain_create("My Domain");

    int num_particles = particles.size();
    Eigen::SparseMatrix<int> interaction_graph(num_particles, num_particles);

    std::vector<double> step_size;
    step_size.resize(num_particles);
    for (int p_i = 0; p_i < num_particles; p_i++) {
        step_size[p_i] = particles[p_i].last_time_step_size;
    }

    std::mutex lock;

    for (int b_i = 0; b_i < broad_phase.size(); b_i++) {
        collision::BroadPhaseCollision& b = broad_phase[b_i];

        int a_id = b.first.first; // geo id
        int b_id = b.second.first;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            // std::lock_guard<std::mutex> guard(lock);
            interaction_graph.coeffRef(a_id, b_id) = 1;
            interaction_graph.coeffRef(b_id, a_id) = 1;
        }
    }

    Eigen::MatrixXi components;
    Eigen::MatrixXi sizes;
    igl::connected_components(interaction_graph, components, sizes);


    cluster_particles_out.clear();
    cluster_interactions_out.clear();
    cluster_step_size_out.clear();
    sleeping.clear();

    std::vector<double> min_current_time;

    for (int c = 0; c < sizes.rows(); c++) {
        cluster_particles_out.emplace_back();
        cluster_particles_out[c].reserve(sizes(c, 0));
        cluster_interactions_out.emplace_back();
        cluster_step_size_out.push_back(max_time);
        sleeping.push_back(true);

        min_current_time.push_back(std::numeric_limits<double>::infinity());
    }

    for (int p_i = 0; p_i < num_particles; p_i++) {
        Particle* address = &(particles[p_i]);
        cluster_particles_out[components(p_i, 0)].push_back(address);
        min_current_time[components(p_i, 0)] = std::min(min_current_time[components(p_i, 0)], particles[p_i].current_state.getTime());
        // cluster_step_size_out[components(p_i, 0)] = std::min(cluster_step_size_out[components(p_i, 0)], step_size[p_i]);
        if (!particles[p_i].getSleeping()) {
            sleeping[components(p_i, 0)] = false;
        }
    }

    for (int c = 0; c < sizes.rows(); c++) {
        for (Particle* p : cluster_particles_out[c]) {
            p->rollBackState(min_current_time[c]);
        }
    }


    // TODO Compute a min future state per cluster given the slow increase of timestep size.  Then use this to avoid computing step sizes for clusters
    //      where the current time is greater than the max possible future time of another cluster.

    #pragma omp parallel
    {
        // for (collision::BroadPhaseCollision& b : broad_phase) {
        #pragma omp single
        {
            #pragma omp taskloop
            for (int b_i = 0; b_i < broad_phase.size(); b_i++) {
                __itt_string_handle* continuous_pair_task = __itt_string_handle_create("Continuous pair");
                __itt_task_begin(domain, __itt_null, __itt_null, continuous_pair_task);

                int tid = omp_get_thread_num();

                collision::BroadPhaseCollision& b = broad_phase[b_i];

                int a_id = b.first.first; // geo id
                int b_id = b.second.first;

                Eigen::Vector3d a_centre = particles[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = particles[b_id].current_state.getTranslation();

                double max_time_step_for_pair = std::min(particles[a_id].last_time_step_size, particles[b_id].last_time_step_size); // Use the minimum of the maximum allowed timesteps

                double min_time;

                __itt_string_handle* continuous_trees_continuous_task = __itt_string_handle_create("Compare trees full continuous");
                __itt_task_begin(domain, __itt_null, __itt_null, continuous_trees_continuous_task);
                std::vector<collision::ContinuousContact<double>> Ccs = collision::compareTreesFullContinuous<double, 8, 8>(particles[a_id], particles[b_id], max_time_step_for_pair, min_time);
                __itt_task_end(domain);

                // if (Ccs.size() > 0) {
                //     step_size[a_id] = std::min(step_size[a_id], min_time * 0.5);
                //     step_size[b_id] = std::min(step_size[b_id], min_time * 0.5);
                // }

                if (Ccs.size() > 0) {
                    double min_scaling_seen = 1.0;
                    double min_toc = 1.0;
                    for (collision::ContinuousContact<double>& c : Ccs) {
                        Eigen::Vector3d hit_normal = (c.A - c.B);

                        double interaction_dist = c.eps_a + c.eps_b;

                        Eigen::Vector<double, 3> a_vel = c.p_a->futurePointVelocity(c.A);
                        Eigen::Vector<double, 3> b_vel = c.p_b->futurePointVelocity(c.B);

                        Eigen::Vector<double, 3> rel_vel = a_vel - b_vel;

                        if (hit_normal.normalized().dot(rel_vel.normalized()) <= 0.0) {
                            if (hit_normal.norm() < 1e-8) {
                                if (c.toc == 0.0) {
                                    throw std::runtime_error("Invalid configuration");
                                    continue;
                                }
                                double rel_vel_norm = rel_vel.norm();
                                double toc_start_contact = std::max(0.0, c.toc - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                // double toc_start_contact = std::max(0.0, 1.0 - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                double new_step = common::lerp(toc_start_contact, c.toc, 0.75);
                                min_scaling_seen = std::min(min_scaling_seen, new_step);
                                // This is only correct if the velocity is perpendicular to the face tangent. Otherwise it's an underestimate of the toc (which is safe).
                                // min_scaling_seen = std::min(min_scaling_seen, c.toc * 0.9);
                                min_toc = std::min(min_toc, c.toc);
                            }
                            else {
                                if (c.toc < 1.0) {
                                    min_toc = std::min(min_toc, c.toc);
                                    // End point is inside range (since it's moving towards each other) => compute start point
                                    double proj_vel = -(hit_normal.normalized().dot(rel_vel)); // +ve
                                    double depth = interaction_dist - hit_normal.norm();
                                    double proj_dist = proj_vel * c.toc * max_time_step_for_pair;
                                    if (proj_dist > 1e-6) {
                                        if (depth < proj_dist) {
                                            // This is the first frame penetrating this eps boundry
                                            assert((2.0 - depth / proj_dist) / 2.0 <= 1.0);

                                            double new_step = common::lerp(c.toc * (1.0 - depth / proj_dist), c.toc, 0.75);
                                            min_scaling_seen = std::min(min_scaling_seen, new_step);
                                            if (max_time_step_for_pair * new_step < 1e-6) {
                                                printf("Small time step 0\n");
                                            }
                                        }
                                        else {
                                            // The contact point was already inside the eps boundry at the start of the timestep
                                            double future_point = c.toc + (1.0 - c.toc) * interaction_dist / (proj_vel * max_time_step_for_pair);

                                            double new_step = common::lerp(c.toc, std::min(future_point, 1.0), 0.9);

                                            min_scaling_seen = std::min(min_scaling_seen, new_step);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    double new_time_step_size = max_time_step_for_pair * min_scaling_seen;
                    double use_time_step_size = std::min(new_time_step_size, std::min(step_size[a_id], step_size[b_id]));

                    std::lock_guard<std::mutex> guard(lock);
                    if (!particles[a_id].is_static) {
                        step_size[a_id] = use_time_step_size;
                    }
                    if (!particles[b_id].is_static) {
                        step_size[b_id] = use_time_step_size;
                    }
                    printf("Min toc: %f\n", min_toc);
                }
                __itt_task_end(domain);
            }
        }
    }

    for (collision::BroadPhaseCollision& b : broad_phase) {
        int a_id = b.first.first; // geo id
        int b_id = b.second.first;
        int id;
        int static_id = -1;
        int cluster_id;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            assert(components(a_id, 0) == components(b_id, 0));
            id = a_id;
        }
        else if (particles[a_id].is_static && particles[b_id].is_static) {
            continue;
        }
        else if (particles[a_id].is_static) {
            id = b_id;
            static_id = a_id;
            cluster_id = components(b_id, 0);
        }
        else {
            id = a_id;
            static_id = b_id;
            cluster_id = components(a_id, 0);
        }

        if (static_id >= 0) {
            Particle* address = &(particles[static_id]);
            if (std::find(cluster_particles_out[cluster_id].begin(), cluster_particles_out[cluster_id].end(), address) == cluster_particles_out[cluster_id].end()) {
                cluster_particles_out[cluster_id].push_back(address);
            }
        }

        cluster_interactions_out[components(a_id, 0)].push_back(b);
        cluster_step_size_out[components(a_id, 0)] = std::min(cluster_step_size_out[components(a_id, 0)], step_size[id]);
    }
}


std::vector<collision::Cluster> collision::separateCollisionClusters(collision::BroadPhaseCollisions& broad_phase, model::ParticleHandler& particles) {
    __itt_domain* domain = __itt_domain_create("My Domain");

    int num_particles = particles.size();
    Eigen::SparseMatrix<int> interaction_graph(num_particles, num_particles);

    for (int b_i = 0; b_i < broad_phase.size(); b_i++) {
        collision::BroadPhaseCollision& b = broad_phase[b_i];

        int a_id = particles.getLocalID(b.first.first); // local geo id
        int b_id = particles.getLocalID(b.second.first);

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            interaction_graph.coeffRef(a_id, b_id) = 1;
            interaction_graph.coeffRef(b_id, a_id) = 1;
        }
    }

    Eigen::MatrixXi components;
    Eigen::MatrixXi sizes;
    igl::connected_components(interaction_graph, components, sizes);

    std::vector<Cluster> clusters;
    std::vector<std::vector<Particle*>> cluster_particles;

    for (int c = 0; c < sizes.rows(); c++) {
        Cluster& C = clusters.emplace_back();

        C.sleeping = true;
        C.min_current_time = std::numeric_limits<double>::infinity();
        C.step_size = std::numeric_limits<double>::infinity();
        C.is_static = true;

        cluster_particles.push_back({});
    }

    for (int p_i = 0; p_i < num_particles; p_i++) {
        Particle* address = &(particles[p_i]);
        cluster_particles[components(p_i, 0)].push_back(address);
        clusters[components(p_i, 0)].min_current_time = std::min(clusters[components(p_i, 0)].min_current_time, particles[p_i].current_state.getTime());
        clusters[components(p_i, 0)].step_size = std::min(clusters[components(p_i, 0)].step_size, address->last_time_step_size);
        if (!particles[p_i].getSleeping()) {
            clusters[components(p_i, 0)].sleeping = false;
        }
    }

    for (int c_i = 0; c_i < sizes.rows(); c_i++) {
        for (Particle* p : cluster_particles[c_i]) {
            p->rollBackState(clusters[c_i].min_current_time);
            if (!p->is_static) {
                clusters[c_i].is_static = false;
            }
        }
    }

    for (collision::BroadPhaseCollision& b : broad_phase) {
        int a_id = particles.getLocalID(b.first.first); // local geo id
        int b_id = particles.getLocalID(b.second.first);
        int id;
        int static_id = -1;
        int cluster_id;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            assert(components(a_id, 0) == components(b_id, 0));
            id = a_id;
        }
        else if (particles[a_id].is_static && particles[b_id].is_static) {
            continue;
        }
        else if (particles[a_id].is_static) {
            id = b_id;
            static_id = a_id;
            cluster_id = components(b_id, 0);
        }
        else {
            id = a_id;
            static_id = b_id;
            cluster_id = components(a_id, 0);
        }

        if (static_id >= 0) {
            Particle* address = &(particles[static_id]);
            if (std::find(cluster_particles[cluster_id].begin(), cluster_particles[cluster_id].end(), address) == cluster_particles[cluster_id].end()) {
                cluster_particles[cluster_id].push_back(address);
            }
        }

        clusters[components(a_id, 0)].interations.push_back(b);
    }

    for (int cluster_i = 0; cluster_i < clusters.size(); cluster_i++) {
        clusters[cluster_i].particles = model::ParticleHandler(cluster_particles[cluster_i]);
    }

    return clusters;
}


void collision::fineCollisionClustersWithTimeStepSelection(Cluster& cluster) {
    int num_particles = cluster.particles.size();
    std::vector<double> step_size;
    step_size.resize(num_particles);
    for (int p_i = 0; p_i < num_particles; p_i++) {
        step_size[p_i] = cluster.particles[p_i].last_time_step_size;
    }

    cluster.step_size = -1.0;

    std::mutex lock;

    #pragma omp parallel
    {
        // for (collision::BroadPhaseCollision& b : broad_phase) {
        #pragma omp single
        {
            #pragma omp taskloop
            for (int b_i = 0; b_i < cluster.interations.size(); b_i++) {
                int tid = omp_get_thread_num();

                collision::BroadPhaseCollision& b = cluster.interations[b_i];

                int a_id = cluster.particles.getLocalID(b.first.first); // geo id
                int b_id = cluster.particles.getLocalID(b.second.first);

                Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                double max_time_step_for_pair = std::min(cluster.particles[a_id].last_time_step_size, cluster.particles[b_id].last_time_step_size); // Use the minimum of the maximum allowed timesteps

                double min_time;

                std::vector<collision::ContinuousContact<double>> Ccs = collision::compareTreesFullContinuous<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], max_time_step_for_pair, min_time);

                if (Ccs.size() > 0) {
                    double min_scaling_seen = 1.0;
                    double min_toc = 1.0;
                    for (collision::ContinuousContact<double>& c : Ccs) {
                        Eigen::Vector3d hit_normal = (c.A - c.B);

                        double interaction_dist = c.eps_a + c.eps_b;

                        Eigen::Vector<double, 3> a_vel = c.p_a->futurePointVelocity(c.A);
                        Eigen::Vector<double, 3> b_vel = c.p_b->futurePointVelocity(c.B);

                        Eigen::Vector<double, 3> rel_vel = a_vel - b_vel;

                        if (hit_normal.normalized().dot(rel_vel.normalized()) <= 0.0) {
                            if (hit_normal.norm() < 1e-8) {
                                if (c.toc == 0.0) {
                                    throw std::runtime_error("Invalid configuration");
                                    continue;
                                }
                                double rel_vel_norm = rel_vel.norm();
                                double toc_start_contact = std::max(0.0, c.toc - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                // double toc_start_contact = std::max(0.0, 1.0 - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                double new_step = common::lerp(toc_start_contact, c.toc, 0.75);
                                min_scaling_seen = std::min(min_scaling_seen, new_step);
                                // This is only correct if the velocity is perpendicular to the face tangent. Otherwise it's an underestimate of the toc (which is safe).
                                // min_scaling_seen = std::min(min_scaling_seen, c.toc * 0.9);
                                min_toc = std::min(min_toc, c.toc);
                            }
                            else {
                                if (c.toc < 1.0) {
                                    min_toc = std::min(min_toc, c.toc);
                                    // End point is inside range (since it's moving towards each other) => compute start point
                                    double proj_vel = -(hit_normal.normalized().dot(rel_vel)); // +ve
                                    double depth = interaction_dist - hit_normal.norm();
                                    double proj_dist = proj_vel * c.toc * max_time_step_for_pair;
                                    if (proj_dist > 1e-6) {
                                        if (depth < proj_dist) {
                                            // This is the first frame penetrating this eps boundry
                                            assert((2.0 - depth / proj_dist) / 2.0 <= 1.0);

                                            double new_step = common::lerp(c.toc * (1.0 - depth / proj_dist), c.toc, 0.75);
                                            min_scaling_seen = std::min(min_scaling_seen, new_step);
                                            if (max_time_step_for_pair * new_step < 1e-6) {
                                                printf("Small time step 0\n");
                                            }
                                        }
                                        else {
                                            // The contact point was already inside the eps boundry at the start of the timestep
                                            double future_point = c.toc + (1.0 - c.toc) * interaction_dist / (proj_vel * max_time_step_for_pair);

                                            double new_step = common::lerp(c.toc, std::min(future_point, 1.0), 0.9);

                                            min_scaling_seen = std::min(min_scaling_seen, new_step);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    double new_time_step_size = max_time_step_for_pair * min_scaling_seen;
                    double use_time_step_size = std::min(new_time_step_size, std::min(step_size[a_id], step_size[b_id]));

                    std::lock_guard<std::mutex> guard(lock);
                    if (!cluster.particles[a_id].is_static) {
                        step_size[a_id] = use_time_step_size;
                    }
                    if (!cluster.particles[b_id].is_static) {
                        step_size[b_id] = use_time_step_size;
                    }
                }
            }
        }
    }

    for (int b_i = 0; b_i < cluster.interations.size(); b_i++) {
        collision::BroadPhaseCollision& b = cluster.interations[b_i];
        int a_id = cluster.particles.getLocalID(b.first.first); // geo id
        int b_id = cluster.particles.getLocalID(b.second.first);
        int id;

        if (!cluster.particles[a_id].is_static && !cluster.particles[b_id].is_static) {
            id = a_id;
        }
        else if (cluster.particles[a_id].is_static && cluster.particles[b_id].is_static) {
            continue;
        }
        else if (cluster.particles[a_id].is_static) {
            id = b_id;
        }
        else {
            id = a_id;
        }

        if (cluster.step_size < 0.0) {
            cluster.step_size = step_size[a_id];
        }
        else {
            cluster.step_size = std::min(cluster.step_size, step_size[id]);
        }
    }
    if (cluster.step_size < 0.0) {
        cluster.step_size = (cluster.particles[0]).last_time_step_size;
    }
}


void collision::fineWitnessCollisionClustersWithTimeStepSelection(std::vector<Delta2::Particle>& particles, std::vector<std::vector<Particle*>>& cluster_particles_out, std::vector<collision::BroadPhaseCollisions>& cluster_interactions_out, std::vector<double>& cluster_step_size_out, std::vector<bool>& sleeping, std::vector<double>& min_current_time, collision::ContactStateCache& cache) {
    __itt_domain* domain = __itt_domain_create("My Domain");

    double max_reachable_time = 0.0;

    for (int c = 0; c < cluster_particles_out.size(); c++) {
        for (Particle* p : cluster_particles_out[c]) {
            max_reachable_time = std::max(max_reachable_time, min_current_time[c] + p->last_time_step_size);
        }
    }

    // Match the current time for all particles in each cluster
    for (int c = 0; c < cluster_particles_out.size(); c++) {
        if (min_current_time[c] < max_reachable_time) {
            for (Particle* p : cluster_particles_out[c]) {
                if (!p->is_static) {
                    p->rollBackState(min_current_time[c]);
                }
            }
        }
    }

    int num_particles = particles.size();
    std::vector<double> step_size;
    step_size.resize(num_particles);
    for (int p_i = 0; p_i < num_particles; p_i++) {
        step_size[p_i] = particles[p_i].last_time_step_size;
    }

    cluster_step_size_out.clear();
    for (int c = 0; c < cluster_particles_out.size(); c++) {
       cluster_step_size_out.push_back(-1.0);
    }

    std::mutex lock;

    #pragma omp parallel
    {
        // for (collision::BroadPhaseCollision& b : broad_phase) {
        #pragma omp single
        {
            #pragma omp taskloop
            for (int c_i = 0; c_i < cluster_interactions_out.size(); c_i++) {
                if (min_current_time[c_i] > max_reachable_time) {
                    continue;
                }
                #pragma omp taskloop
                for (int b_i = 0; b_i < cluster_interactions_out[c_i].size(); b_i++) {
                    __itt_string_handle* continuous_pair_task = __itt_string_handle_create("Continuous pair");
                    __itt_task_begin(domain, __itt_null, __itt_null, continuous_pair_task);

                    int tid = omp_get_thread_num();

                    collision::BroadPhaseCollision& b = cluster_interactions_out[c_i][b_i];

                    int a_id = b.first.first; // geo id
                    int b_id = b.second.first;

                    Eigen::Vector3d a_centre = particles[a_id].current_state.getTranslation();
                    Eigen::Vector3d b_centre = particles[b_id].current_state.getTranslation();

                    double max_time_step_for_pair = std::min(particles[a_id].last_time_step_size, particles[b_id].last_time_step_size); // Use the minimum of the maximum allowed timesteps

                    double min_time;

                    __itt_string_handle* continuous_trees_continuous_task = __itt_string_handle_create("Compare trees full continuous");
                    __itt_task_begin(domain, __itt_null, __itt_null, continuous_trees_continuous_task);
                    // std::vector<collision::ContinuousContact<double>> Ccs = collision::compareTreesFullContinuous<double, 8, 8>(particles[a_id], particles[b_id], max_time_step_for_pair, min_time);
                    std::vector<collision::ContinuousContact<double>> Ccs = collision::compareTreesFullContinuousWitnesses(particles[a_id], particles[b_id], cache, max_time_step_for_pair, min_time);

                    __itt_task_end(domain);

                    // if (Ccs.size() > 0) {
                    //     step_size[a_id] = std::min(step_size[a_id], min_time * 0.5);
                    //     step_size[b_id] = std::min(step_size[b_id], min_time * 0.5);
                    // }

                    if (Ccs.size() > 0) {
                        double min_scaling_seen = 1.0;
                        double min_toc = 1.0;
                        for (collision::ContinuousContact<double>& c : Ccs) {
                            Eigen::Vector3d hit_normal = (c.A - c.B);

                            double interaction_dist = c.eps_a + c.eps_b;

                            Eigen::Vector<double, 3> a_vel = c.p_a->futurePointVelocity(c.A);
                            Eigen::Vector<double, 3> b_vel = c.p_b->futurePointVelocity(c.B);

                            Eigen::Vector<double, 3> rel_vel = a_vel - b_vel;

                            if (hit_normal.normalized().dot(rel_vel.normalized()) <= 0.0) {
                                if (hit_normal.norm() < 1e-8) {
                                    if (c.toc == 0.0) {
                                        throw std::runtime_error("Invalid configuration");
                                        continue;
                                    }
                                    double rel_vel_norm = rel_vel.norm();
                                    double toc_start_contact = std::max(0.0, c.toc - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                    // double toc_start_contact = std::max(0.0, 1.0 - (rel_vel.norm() * max_time_step_for_pair) / interaction_dist);
                                    double new_step = common::lerp(toc_start_contact, c.toc, 0.75);
                                    min_scaling_seen = std::min(min_scaling_seen, new_step);
                                    // This is only correct if the velocity is perpendicular to the face tangent. Otherwise it's an underestimate of the toc (which is safe).
                                    // min_scaling_seen = std::min(min_scaling_seen, c.toc * 0.9);
                                    min_toc = std::min(min_toc, c.toc);
                                }
                                else {
                                    if (c.toc < 1.0) {
                                        min_toc = std::min(min_toc, c.toc);
                                        // End point is inside range (since it's moving towards each other) => compute start point
                                        double proj_vel = -(hit_normal.normalized().dot(rel_vel)); // +ve
                                        double depth = interaction_dist - hit_normal.norm();
                                        double proj_dist = proj_vel * c.toc * max_time_step_for_pair;
                                        if (proj_dist > 1e-6) {
                                            if (depth < proj_dist) {
                                                // This is the first frame penetrating this eps boundry
                                                assert((2.0 - depth / proj_dist) / 2.0 <= 1.0);

                                                double new_step = common::lerp(c.toc * (1.0 - depth / proj_dist), c.toc, 0.75);
                                                min_scaling_seen = std::min(min_scaling_seen, new_step);
                                                if (max_time_step_for_pair * new_step < 1e-6) {
                                                    printf("Small time step 0\n");
                                                }
                                            }
                                            else {
                                                // The contact point was already inside the eps boundry at the start of the timestep
                                                double future_point = c.toc + (1.0 - c.toc) * interaction_dist / (proj_vel * max_time_step_for_pair);

                                                double new_step = common::lerp(c.toc, std::min(future_point, 1.0), 0.9);

                                                min_scaling_seen = std::min(min_scaling_seen, new_step);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        double new_time_step_size = max_time_step_for_pair * min_scaling_seen;
                        double use_time_step_size = std::min(new_time_step_size, std::min(step_size[a_id], step_size[b_id]));

                        std::lock_guard<std::mutex> guard(lock);
                        if (!particles[a_id].is_static) {
                            step_size[a_id] = use_time_step_size;
                        }
                        if (!particles[b_id].is_static) {
                            step_size[b_id] = use_time_step_size;
                        }
                    }
                    __itt_task_end(domain);
                }
            }
        }
    }

    for (int c_i = 0; c_i < cluster_interactions_out.size(); c_i++) {
        for (int b_i = 0; b_i < cluster_interactions_out[c_i].size(); b_i++) {
            collision::BroadPhaseCollision& b = cluster_interactions_out[c_i][b_i];
            int a_id = b.first.first; // geo id
            int b_id = b.second.first;
            int id;

            if (!particles[a_id].is_static && !particles[b_id].is_static) {
                id = a_id;
            }
            else if (particles[a_id].is_static && particles[b_id].is_static) {
                continue;
            }
            else if (particles[a_id].is_static) {
                id = b_id;
            }
            else {
                id = a_id;
            }

            if (cluster_step_size_out[c_i] < 0.0) {
                cluster_step_size_out[c_i] = step_size[a_id];
            }
            else {
                cluster_step_size_out[c_i] = std::min(cluster_step_size_out[c_i], step_size[id]);
            }
        }
    }
    for (int c_i = 0; c_i < cluster_step_size_out.size(); c_i++) {
        if (cluster_step_size_out[c_i] < 0.0) {
            cluster_step_size_out[c_i] = cluster_particles_out[c_i][0]->last_time_step_size;
        }
    }
}


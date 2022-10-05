#include "explicit_adaptive_clusters.h"
#include "../collision_detection/broad_phase_embree.h"
#include "../collision_detection/full_tree_comparison.h"
#include "../model/forces.h"
#include "../model/integrate.h"
#include "../collision_detection/separate_clusters.h"

#include <chrono>
#include <omp.h>
#include <ittnotify.h>

using namespace Delta2;

struct ExplicitState
{
    ExplicitState()
    {
        force = {0, 0, 0};
        torque = {0, 0, 0};
    }
    Eigen::Vector3d force;
    Eigen::Vector3d torque;
};

timestepping::ExplicitAdaptiveClustersScheme::ExplicitAdaptiveClustersScheme(std::vector<Particle> *ps, std::function<Eigen::Vector3d(const Particle &)> external_force, common::Options& o) : timestepping::TimestepScheme(ps, external_force, o)
{
    _last_time_step_size = -1.0;

    int id = 1;

    _sleep_candidates.resize(_particles->size());
    for (int i = 0; i < _particles->size(); i++) {
        _sleep_candidates[i] = 0.0;
    }

    for (Particle &p : *_particles)
    {
        p.last_state = p.current_state;
    }
};

void timestepping::ExplicitAdaptiveClustersScheme::step(double time, std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> &view_draws)
{
	__itt_domain* domain = __itt_domain_create("My Domain");

    auto t0 = std::chrono::high_resolution_clock::now();
    __itt_string_handle* setup_task = __itt_string_handle_create("Setup particles");
    __itt_task_begin(domain, __itt_null, __itt_null, setup_task);

    if (_last_time_step_size <= 0.0)
    {
        _last_time_step_size = time;

        for (Particle &p : *_particles)
        {
            p.last_time_step_size = time;
        }
    }

    for (Particle &p : *_particles)
    {
        p.last_time_step_size = std::min(p.last_time_step_size * 1.1, time);
        p.projectFutureState(p.last_time_step_size);
    }

    __itt_task_end(domain);
    auto t1 = std::chrono::high_resolution_clock::now();
    printf("Setup particles:              %ldms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    __itt_string_handle* broad_phase_task = __itt_string_handle_create("Broad phase");
    __itt_task_begin(domain, __itt_null, __itt_null, broad_phase_task);

    collision::BroadPhaseCollisions B;
    B = collision::broadPhaseEmbree(*_particles);

    __itt_task_end(domain);
    auto t2 = std::chrono::high_resolution_clock::now();
    printf("Broad phase:                  %ldms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    __itt_string_handle* separate_task = __itt_string_handle_create("Separate");
    __itt_task_begin(domain, __itt_null, __itt_null, separate_task);

    std::vector<double> min_times;
    double min_final_time = std::numeric_limits<double>::max();

    std::vector<std::vector<Delta2::Particle *>> cluster_particles;
    std::vector<collision::BroadPhaseCollisions> cluster_interactions;
    std::vector<double> cluster_step_size;
    std::vector<bool> cluster_sleeping;
    // collision::separateCollisionClustersWithTimeStepSelection(B, *_particles, time, cluster_particles, cluster_interactions, cluster_step_size, cluster_sleeping);
    
    std::vector<double> min_current_times;
    collision::separateCollisionClusters(B, *_particles, cluster_particles, cluster_interactions, cluster_sleeping, min_current_times);
    collision::fineCollisionClustersWithTimeStepSelection(*_particles, cluster_particles, cluster_interactions, cluster_step_size, cluster_sleeping, min_current_times);
    
    __itt_task_end(domain);
    auto t2_5 = std::chrono::high_resolution_clock::now();
    printf("Separate Clusters:            %ldms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2_5 - t2).count());
    __itt_string_handle* min_times_task = __itt_string_handle_create("Min times");
    __itt_task_begin(domain, __itt_null, __itt_null, min_times_task);
    
    for (int cluster_i = 0; cluster_i < cluster_particles.size(); cluster_i++)
    {
        // Prevent the time step size going to zero even if there is going to be an intersection
        cluster_step_size[cluster_i] = std::max(cluster_step_size[cluster_i], double(_opt.time_step_size * _opt.adaptive_time_step_factor));

        double min_time = std::numeric_limits<double>::max();
        bool has_non_static = false;
        for (Particle *p : cluster_particles[cluster_i])
        {
            if (!p->is_static)
            {
                min_time = std::min(min_time, p->current_state.getTime());
                has_non_static = true;
            }
        }
        if (has_non_static)
        {
            min_times.emplace_back(min_time);
            printf("cluster %i min time: %f\n", cluster_i, min_time);
            // for (Particle *p : cluster_particles[cluster_i])
            // {
            //     printf("    Particle %i\n", p->id);
            // }
            min_final_time = std::min(min_final_time, min_time + cluster_step_size[cluster_i]);
        }
        else
        {
            min_times.emplace_back(-1);
        }
    }

    __itt_task_end(domain);
    auto t3 = std::chrono::high_resolution_clock::now();
    printf("particlesMin times:                    %ldms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2_5).count());
    __itt_string_handle* process_clusters_task = __itt_string_handle_create("Process all clusters");
    __itt_task_begin(domain, __itt_null, __itt_null, process_clusters_task);
    __itt_string_handle* process_cluster_task = __itt_string_handle_create("Process cluster");

    for (int cluster_i = 0; cluster_i < cluster_particles.size(); cluster_i++)
    {
        __itt_task_begin(domain, __itt_null, __itt_null, process_cluster_task);

        if (min_times[cluster_i] > min_final_time)
        {
            __itt_task_end(domain);
            continue;
        }

        if (cluster_sleeping[cluster_i] || min_times[cluster_i] == -1) {
            __itt_task_end(domain);
            continue;
        }

        std::vector<Eigen::Vector3d> forces;
        std::vector<Eigen::Vector3d> torques;
        std::vector<int> counts;
        std::vector<bool> can_advance;

        std::unordered_map<int, int> global_to_local_ids;

        for (int p_i = 0; p_i < cluster_particles[cluster_i].size(); p_i++)
        {
            Particle *p = cluster_particles[cluster_i][p_i];
            forces.emplace_back(0.0, 0.0, 0.0);
            torques.emplace_back(0.0, 0.0, 0.0);
            counts.emplace_back(0);

            global_to_local_ids[p->id] = p_i;
            // printf("%i: %i\n", p->id, p_i);
        }

        __itt_string_handle* compare_trees_task = __itt_string_handle_create("Compare all trees");
        __itt_task_begin(domain, __itt_null, __itt_null, compare_trees_task);
        std::vector<std::tuple<int, int, collision::Contact<double>>> hits;
        bool valid_timestep_found = false;
        while (!valid_timestep_found)
        {
            printf("cluster %i time step size: %f\n", cluster_i, cluster_step_size[cluster_i]);
            valid_timestep_found = true;
            hits.clear();

            for (Particle *p : cluster_particles[cluster_i])
            {
                p->last_time_step_size = cluster_step_size[cluster_i];
                if (!p->is_static)
                {
                    p->rollBackState(min_times[cluster_i]);
                    p->projectFutureState(p->last_time_step_size);
                }
                else
                {
                    p->current_state.setTime(min_times[cluster_i]);
                    p->future_state = p->current_state;
                    p->future_state.setTime(p->last_time_step_size);
                }
            }

            std::mutex lock;

            //for (collision::BroadPhaseCollision &b : cluster_interactions[cluster_i])
            // #pragma omp parallel
            {
                // #pragma omp single
                {
                    // #pragma omp taskloop
                    for (int b_i = 0; b_i < cluster_interactions[cluster_i].size(); b_i++)
                    {
                        __itt_string_handle* tree_comparison_task = __itt_string_handle_create("Tree comparison");
                        __itt_task_begin(domain, __itt_null, __itt_null, tree_comparison_task);
                        int tid = omp_get_thread_num();

                        collision::BroadPhaseCollision &b = cluster_interactions[cluster_i][b_i]; 

                        int a_id = b.first.first; // geo id
                        int b_id = b.second.first;

                        Eigen::Vector3d a_centre = (*_particles)[a_id].current_state.getTranslation();
                        Eigen::Vector3d b_centre = (*_particles)[b_id].current_state.getTranslation();

                        std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>((*_particles)[a_id], (*_particles)[b_id]);
                        std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

                        std::lock_guard<std::mutex> guard(lock);
                        for (collision::Contact<double> &c : filtered)
                        {
                            if ((c.A - c.B).norm() < 1e-8)
                            {
                                if (cluster_step_size[cluster_i] * 0.9 > double(_opt.time_step_size * _opt.adaptive_time_step_factor)) {
                                    valid_timestep_found = false;
                                    cluster_step_size[cluster_i] *= 0.9;
                                    break;
                                }
                            }
                            hits.push_back(std::make_tuple(a_id, b_id, c));
                        }
                        __itt_task_end(domain);
                    }
                }
            }
        }
        __itt_task_end(domain);

        __itt_string_handle* accumulate_hits_task = __itt_string_handle_create("Accumulate hits");
        __itt_task_begin(domain, __itt_null, __itt_null, accumulate_hits_task);
        for (auto &[a_id, b_id, c] : hits)
        {
            Eigen::Vector3d a_centre = (*_particles)[a_id].current_state.getTranslation();
            Eigen::Vector3d b_centre = (*_particles)[b_id].current_state.getTranslation();

            Eigen::Vector3d a_hit_point = c.A.cast<double>();
            Eigen::Vector3d b_hit_point = c.B.cast<double>();

            double mass_multiply = 0.0;
            if (!(*_particles)[a_id].is_static)
            {
                mass_multiply += (*_particles)[a_id].getMass();
            }
            if ((*_particles)[b_id].is_static)
            {
                mass_multiply += (*_particles)[b_id].getMass();
            }
            // Eigen::Vector3d f = model::calcForce(c, 1.0, (*_particles)[a_id].last_time_step_size).cast<double>() * mass_multiply;
            Eigen::Vector3d f = model::calcImpulse(c, (*_particles)[a_id].last_time_step_size).cast<double>() / time;

            if (f.norm() > 1e-8)
            {
                Eigen::Vector3d t_a = model::calcTorque(f, a_hit_point, a_centre);
                Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), b_hit_point, b_centre);

                // printf("Force: %f\n", f.norm());

                forces[global_to_local_ids[a_id]] += f;
                torques[global_to_local_ids[a_id]] += t_a;
                counts[global_to_local_ids[a_id]]++;
                forces[global_to_local_ids[b_id]] += -f;
                torques[global_to_local_ids[b_id]] += t_b;
                counts[global_to_local_ids[b_id]]++;

                view_draws.push_back(std::make_pair(a_hit_point, b_hit_point));
            }
        }
        __itt_task_end(domain);
        
        __itt_string_handle* integrate_1st_task = __itt_string_handle_create("Integrate 1st");
        __itt_task_begin(domain, __itt_null, __itt_null, integrate_1st_task);
        for (Particle *p : cluster_particles[cluster_i])
        {
            if (!p->is_static)
            {
                Eigen::Vector3d force = {0.0, 0.0, 0.0};
                Eigen::Vector3d torque = {0.0, 0.0, 0.0};
                if (counts[global_to_local_ids[p->id]])
                {
                    force = forces[global_to_local_ids[p->id]] / double(counts[global_to_local_ids[p->id]]);
                    torque = torques[global_to_local_ids[p->id]] / double(counts[global_to_local_ids[p->id]]);
                }
                model::integrateEuler(*p, force + _external_force(*p), torque, cluster_step_size[cluster_i], false);
            }
            else
            {
                // model::staticApply(*p, cluster_step_size[cluster_i], false);
            }
        }
        __itt_task_end(domain);

        __itt_string_handle* friction_solve_task = __itt_string_handle_create("Friction solve");
        __itt_task_begin(domain, __itt_null, __itt_null, friction_solve_task);
        Delta2::model::friction_solve(*_particles, hits, forces, torques, counts, _external_force, cluster_step_size[cluster_i], global_to_local_ids);
        __itt_task_end(domain);

        __itt_string_handle* integrate_2nd_task = __itt_string_handle_create("Integrate 2nd");
        __itt_task_begin(domain, __itt_null, __itt_null, integrate_2nd_task);
        for (Particle *p : cluster_particles[cluster_i])
        {
            if (!p->is_static)
            {
                Eigen::Vector3d force = {0.0, 0.0, 0.0};
                Eigen::Vector3d torque = {0.0, 0.0, 0.0};
                if (counts[global_to_local_ids[p->id]])
                {
                    force = forces[global_to_local_ids[p->id]] / double(counts[global_to_local_ids[p->id]]);
                    torque = torques[global_to_local_ids[p->id]] / double(counts[global_to_local_ids[p->id]]);
                }
                double diff = model::integrateEuler(*p, force + _external_force(*p), torque, cluster_step_size[cluster_i]);
                if (diff < 2.0 * p->geo_eps) {
                    if (p->current_state.getTime() - _sleep_candidates[p->id] > 0.25) {
                        p->setSleeping(true);
                        printf("Particle: %i put to sleep\n", p->id);
                    }
                }
                else {
                    if (p->getSleeping()) {
                        p->setSleeping(false);
                        printf("Particle: %i woken up\n", p->id);
                    }
                    _sleep_candidates[p->id] = p->current_state.getTime();
                }
            }
            else
            {
                // model::staticApply(*p, cluster_step_size[cluster_i]);
            }
            printf("Post integration state of %i: last time %f, time %f, future time %f, sleep from %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), _sleep_candidates[p->id]);
        }
        __itt_task_end(domain);
        __itt_task_end(domain);
    }

    __itt_task_end(domain);
    auto t4 = std::chrono::high_resolution_clock::now();
    printf("Process Clusters:             %ldms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count());
}

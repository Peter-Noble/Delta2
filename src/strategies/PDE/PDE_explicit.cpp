#include "PDE_explicit.h"

#include "../../collision_detection/full_tree_comparison.h"
#include "../../model/forces.h"

#include <mutex>

using namespace Delta2;
using namespace strategy;

PDEExplicit::PDEExplicit(ContactDetectionStrategy& contact_detection,
                         ContactForceStrategy& contact_force,
                         FrictionStrategy& friction,
                         TimeStepSelectionStrategy& time_step,
                         common::Options& opt) :
                         PDEStrategy(opt),
                         _contact_detection(contact_detection),
                         _contact_force(contact_force),
                         _friction(friction),
                         _time_step(time_step) {
}

double PDEExplicit::selectTimeStep(collision::Cluster& cluster) {
    // printf("Selecting explicit time step\n");
    double time = _time_step.selectTimeStep(cluster);
    cluster.step_size = time;
    return time;
}

bool PDEExplicit::step(collision::Cluster& cluster) {
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<int> counts;
    std::vector<bool> can_advance;

    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    for (Particle* p : cluster.particles) {
        forces.emplace_back(0.0, 0.0, 0.0);
        torques.emplace_back(0.0, 0.0, 0.0);
        counts.emplace_back(0);
    }

    std::vector<collision::Contact<double>> hits;

    for (Particle *p : cluster.particles)
    {
        p->last_time_step_size = cluster.step_size;
        if (!p->is_static)
        {
            p->rollBackState(cluster.min_current_time);
            p->projectFutureState(p->last_time_step_size);
        }
        // else
        // {
        //     p->current_state.setTime(min_times[cluster_i]);
        //     p->future_state = p->current_state;
        //     p->future_state.setTime(p->last_time_step_size);
        // }
    }

    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    std::mutex lock;

    // #pragma omp parallel
    {
        // #pragma omp single
        {
            // #pragma omp taskloop
            for (int b_i = 0; b_i < cluster.interations.size(); b_i++)
            {
                collision::BroadPhaseCollision &b = cluster.interations[b_i]; 

                int a_id = cluster.particles.getLocalID(b.first.first); // geo id
                int b_id = cluster.particles.getLocalID(b.second.first);

                Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
                Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

                std::vector<collision::Contact<double>> Cs = collision::compareTreesFull<double, 8, 8>(cluster.particles[a_id], cluster.particles[b_id], _contact_detection);
                std::vector<collision::Contact<double>> filtered = collision::filterContacts<double>(Cs, 0.0);

                std::lock_guard<std::mutex> guard(lock);
                for (collision::Contact<double> &c : filtered)
                {
                    hits.push_back(c);
                }
            }
        }
    }
    
    bool success = _contact_force.solve(cluster, hits);
    if (!success) {
        return false;
    }
    
    for (Particle* p : cluster.particles)
    {
        assert(p->last_state.getTime() <= p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    bool willSleep = true;
    for (Particle *p : cluster.particles)
    {
        if (!p->is_static)
        {
            // assert(p->id < _not_moved_for.size());
            bool canSleep = p->last_state.isStationary() && p->current_state.isStationary() && p->future_state.isStationary();
            canSleep &= p->last_state.getTime() < p->current_state.getTime() && p->current_state.getTime() < p->future_state.getTime();
            if (!canSleep) {
                p->sleep_not_moved_steps = 0;
            } else {
                p->sleep_not_moved_steps++;
            }
            if (p->sleep_not_moved_steps < 20) {
                willSleep = false;
            }

            // printf("Updating %i last: %f, current: %f, future: %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime());
            p->last_state = p->current_state;
            p->current_state = p->future_state;
            p->projectFutureState(cluster.step_size);

            // printf("Post integration state of %i: last time %f, time %f, future time %f, sleep from %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), p->sleep_candidate_time);
        }
    }

    for (Particle* p : cluster.particles)
    {
        p->setSleeping(willSleep);
        assert(p->last_state.getTime() < p->current_state.getTime() || p->current_state.getTime() == 0 || p->is_static);
    }

    return true;
}

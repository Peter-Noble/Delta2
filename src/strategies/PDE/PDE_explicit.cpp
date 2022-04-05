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
    double time = _time_step.selectTimeStep(cluster);
    cluster.step_size = time;
    return time;
}

void PDEExplicit::step(collision::Cluster& cluster) {
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<int> counts;
    std::vector<bool> can_advance;

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

    std::mutex lock;

    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp taskloop
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
    
    _contact_force.solve(cluster, hits);
    
    for (Particle *p : cluster.particles)
    {
        if (!p->is_static)
        {
            p->last_state = p->current_state;
            p->current_state = p->future_state;
            p->projectFutureState(cluster.step_size);

            printf("Post integration state of %i: last time %f, time %f, future time %f, sleep from %f\n", p->id, p->last_state.getTime(), p->current_state.getTime(), p->future_state.getTime(), p->sleep_candidate_time);
        }
    }
}

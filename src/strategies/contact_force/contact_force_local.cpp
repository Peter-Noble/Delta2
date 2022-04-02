#include "contact_force_local.h"

#include "../../model/forces.h"
#include "../../model/integrate.h"

using namespace Delta2;
using namespace strategy;

ContactForceLocal::ContactForceLocal(FrictionStrategy& friction,
                                     common::Options& opt) :
                                     ContactForceStrategy(friction, opt) {

}

void ContactForceLocal::solve(collision::Cluster& cluster, std::vector<collision::Contact<double>>& hits) {
    std::vector<Eigen::Vector3d> forces;
    std::vector<Eigen::Vector3d> torques;
    std::vector<int> counts;

    for (int p_i = 0; p_i < cluster.particles.size(); p_i++) {
        forces.push_back({0, 0, 0});
        torques.push_back({0, 0, 0});
        counts.push_back(0);
    }
    
    for (collision::Contact<double>& c : hits)
    {
        int a_id = cluster.particles.getLocalID(c.p_a->id);
        int b_id = cluster.particles.getLocalID(c.p_b->id);

        Eigen::Vector3d a_centre = cluster.particles[a_id].current_state.getTranslation();
        Eigen::Vector3d b_centre = cluster.particles[b_id].current_state.getTranslation();

        Eigen::Vector3d a_hit_point = c.A.cast<double>();
        Eigen::Vector3d b_hit_point = c.B.cast<double>();

        double mass_multiply = 0.0;
        if (!cluster.particles[a_id].is_static)
        {
            mass_multiply += cluster.particles[a_id].getMass();
        }
        if (cluster.particles[b_id].is_static)
        {
            mass_multiply += cluster.particles[b_id].getMass();
        }
        // Eigen::Vector3d f = model::calcForce(c, 1.0, cluster.particles[a_id].last_time_step_size).cast<double>() * mass_multiply;
        Eigen::Vector3d f = model::calcImpulse(c, cluster.particles[a_id].last_time_step_size).cast<double>() / cluster.step_size;

        if (f.norm() > 1e-8)
        {
            Eigen::Vector3d t_a = model::calcTorque(f, a_hit_point, a_centre);
            Eigen::Vector3d t_b = model::calcTorque(Eigen::Vector3d(-f), b_hit_point, b_centre);

            // printf("Force: %f\n", f.norm());

            forces[a_id] += f;
            torques[a_id] += t_a;
            counts[a_id]++;
            forces[b_id] += -f;
            torques[b_id] += t_b;
            counts[b_id]++;

            // view_draws.push_back(std::make_pair(a_hit_point, b_hit_point));
        }
    }
    
    for (Particle *p : cluster.particles)
    {
        if (!p->is_static)
        {
            Eigen::Vector3d force = {0.0, 0.0, 0.0};
            Eigen::Vector3d torque = {0.0, 0.0, 0.0};
            uint32_t id = cluster.particles.getLocalID(p->id);
            if (counts[id] > 0)
            {
                force = forces[id] / double(counts[id]);
                torque = torques[id] / double(counts[id]);
            }
            model::integrateEuler(*p, force + external_force(*p), torque, cluster.step_size, false);
        }
    }
    _friction.solve(cluster.particles, hits, forces, torques, counts, external_force, cluster.step_size);
    // Delta2::model::friction_solve(cluster.particles, hits, forces, torques, counts, external_force, cluster.step_size);

    for (Particle *p : cluster.particles)
    {
        if (!p->is_static)
        {
            Eigen::Vector3d force = {0.0, 0.0, 0.0};
            Eigen::Vector3d torque = {0.0, 0.0, 0.0};
            uint32_t id = cluster.particles.getLocalID(p->id);
            if (counts[id] > 0)
            {
                force = forces[id] / double(counts[id]);
                torque = torques[id] / double(counts[id]);
            }
            double diff = model::integrateEuler(*p, force + external_force(*p), torque, cluster.step_size, false);
            if (diff < 2.0 * p->geo_eps) {
                if (p->current_state.getTime() - p->sleep_candidate_time > 0.25) {
                    p->setSleeping(true);
                    printf("Particle: %i put to sleep\n", p->id);
                }
            }
            else {
                if (p->getSleeping()) {
                    p->setSleeping(false);
                    printf("Particle: %i woken up\n", p->id);
                }
                p->sleep_candidate_time = p->current_state.getTime();
            }
        }
    }
}

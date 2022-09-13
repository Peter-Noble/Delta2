#pragma once

#include "../model/particle_handler.h"
#include "broad_phase.h"


namespace Delta2 {
    namespace collision {
        struct Cluster {
            model::ParticleHandler particles;
            collision::BroadPhaseCollisions interations;
            double step_size;
            bool sleeping;
            double min_current_time;
            bool is_static;

            // Check to see if this other cluster is a match for this one and hasn't advanced in time.
            bool hasAdvanced(const Cluster& other) const {
                if (min_current_time != other.min_current_time) {
                    return true;
                }
                if (sleeping != other.sleeping) {
                    return true;
                }
                if (particles.size() != other.particles.size()) {
                    return true;
                }
                for (Particle* p : particles) {
                    bool found = false;
                    for (const Particle* o : other.particles) {
                        if (p->id == o->id) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        return true;
                    }
                }
            };
        };
    }
}
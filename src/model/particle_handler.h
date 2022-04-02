#pragma once

#include "particle.h"

namespace Delta2 {
    namespace model {
        class ParticleHandler {
        public:
            ParticleHandler(); // TODO Here for convenience so that clusters can be constructed without particles ready yet.  Is there a better way of doingn this?
            ParticleHandler(std::vector<Particle*>& ps);
            ParticleHandler(std::vector<Particle>& ps);
            Particle& getParticle(uint32_t local_id) const;
            Particle& operator[](uint32_t local_id) const;
            Particle& getParticleGlobalID(uint32_t id) const;
            uint32_t getLocalID(uint32_t id) const;
            uint32_t size() const;
            std::vector<Particle*>::iterator begin();
            std::vector<Particle*>::iterator end();
        private:
            std::vector<Particle*> _ps;
            std::unordered_map<uint32_t, uint32_t> _global_to_local_ids;
        };
    }
}
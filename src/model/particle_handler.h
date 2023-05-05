#pragma once

#include "particle.h"

namespace Delta2 {
    namespace model {
        class ParticleHandler {
        public:
            ParticleHandler(); // TODO Here for convenience so that clusters can be constructed without particles ready yet.  Is there a better way of doingn this?
            ParticleHandler(std::vector<Particle*>& ps);
            ParticleHandler(std::vector<Particle>& ps);
            inline Particle& getParticle(uint32_t local_id) const {return *_ps[local_id];};
            inline Particle& operator[](uint32_t local_id) const {return getParticle(local_id);};
            Particle& getParticleGlobalID(uint32_t id) const;
            inline uint32_t getLocalID(uint32_t id) const {return _global_to_local_ids.at(id);};
            uint32_t size() const;
            std::vector<Particle*>::iterator begin();
            std::vector<Particle*>::iterator end();
            std::vector<Particle*>::const_iterator begin() const;
            std::vector<Particle*>::const_iterator end() const;
            bool isValid() const;
            void initLast();
        private:
            std::vector<Particle*> _ps;
            std::unordered_map<uint32_t, uint32_t> _global_to_local_ids;
        };
    }
}
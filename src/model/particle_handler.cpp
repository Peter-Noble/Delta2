#include "particle_handler.h"

using namespace Delta2;
using namespace model;

ParticleHandler::ParticleHandler() {
    
}

ParticleHandler::ParticleHandler(std::vector<Particle*>& ps) {
    _ps = ps;

    for (int p_i = 0; p_i < ps.size(); p_i++) {
        Particle *p = ps[p_i];
        _global_to_local_ids[p->id] = p_i;
    }

    if (!isValid()) {
        throw std::runtime_error("Particle handler doesn't have a valid map!");
    }
}

ParticleHandler::ParticleHandler(std::vector<Particle>& ps) {
    _ps.reserve(ps.size());

    for (int p_i = 0; p_i < ps.size(); p_i++) {
        _ps.push_back(&ps[p_i]);
        _global_to_local_ids[ps[p_i].id] = p_i;
    }

    if (!isValid()) {
        throw std::runtime_error("Particle handler doesn't have a valid map!");
    }
}

Particle& ParticleHandler::getParticle(uint32_t local_id) const {
    return *_ps[local_id];
}

Particle& ParticleHandler::operator[](uint32_t local_id) const {
    return getParticle(local_id);
}

Particle& ParticleHandler::getParticleGlobalID(uint32_t id) const {
    uint32_t local_id = _global_to_local_ids.at(id);
    return *_ps[local_id];
}

uint32_t ParticleHandler::getLocalID(uint32_t id) const {
    return _global_to_local_ids.at(id);
}

uint32_t ParticleHandler::size() const {
    return _ps.size();
}

std::vector<Particle*>::iterator ParticleHandler::begin() {
    return _ps.begin();
}

std::vector<Particle*>::iterator ParticleHandler::end() {
    return _ps.end();
}

std::vector<Particle*>::const_iterator ParticleHandler::begin() const {
    return _ps.cbegin();
}

std::vector<Particle*>::const_iterator ParticleHandler::end() const {
    return _ps.cend();
}

bool ParticleHandler::isValid() const {
    bool result = true;
    for (Particle* p : _ps) {
        if (_global_to_local_ids.find(p->id) == _global_to_local_ids.end()) {
            result = false;
            break;
        }
    }
    assert(result);
    return result;
}

void ParticleHandler::initLast() {
    for (Particle* p : _ps) {
        p->last_state = p->current_state;
    }
}

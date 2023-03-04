#pragma once

#include <vector>
#include <Eigen/Dense>
#include "../model/surrogate_components.h"
#include "../model/particle.h"
#include "../model/collision_state.h"

// #include <bits/stdc++.h>
  

namespace Delta2 {
    namespace collision {
        class SeparationWitness {
        public:
            SeparationWitness(uint32_t a_n, uint32_t b_n, const Eigen::Vector3d& origin, const Eigen::Vector3d& normal, bool a_def);
            SeparationWitness(Particle& a, Particle& b, const DeferredCompare& dc, bool& success);
            bool check(Particle& a, Particle& b);
            DeferredCompare toDeferredCompare(Particle& a, Particle& b);
        private:
            Eigen::Vector3d origin;
            Eigen::Vector3d normal; // Normal points away from the defining
            uint32_t a_node;
            uint32_t b_node;
            bool a_defining;
            double eps;
        };

        class ContactState {
        public:
            ContactState();
            ContactState(DeferredCompare root);
            bool attemptWitness(Particle& a, Particle& b, const DeferredCompare& dc);
            std::vector<SeparationWitness> witnesses;
            std::vector<DeferredCompare> noWitnesses;
        };

        // A hash function used to hash a pair of any kind
        struct HashContactPair {
            template <class T1, class T2>
            size_t operator()(const std::pair<T1, T2>& p) const
            {
                auto hash1 = std::hash<T1>{}(p.first);
                auto hash2 = std::hash<T2>{}(p.second);
                return hash1 ^ hash2;
            }
        };

        class ContactStateCache {
        public:
            ContactStateCache();
            ContactState& getContactState(Particle& a_id, Particle& b_id);
        private:
            std::unordered_map<std::pair<uint32_t, uint32_t>, ContactState, HashContactPair> states;
            std::mutex _lock;
        };
    }
}

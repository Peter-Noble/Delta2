#pragma once

#include "LocalContact.h"

#include <unordered_map>
#include <mutex>
#include <shared_mutex>

namespace Delta2 {
    struct index_pair_hash {
        std::size_t operator () (std::pair<int, int> const &v) const
        {
            return std::hash<int>()(v.first + v.second * 100000);    
        }
    };
    class SequentialImpulseWarmStart {
        public:
            SequentialImpulseWarmStart();
            void apply(std::vector<LocalContact>& contacts, std::vector<collision::Contact<double>>& hits);
            void update(std::vector<LocalContact>& contacts, std::vector<collision::Contact<double>>& hits);
        private:
            std::unordered_map<std::pair<int, int>, double, index_pair_hash> _store;
            std::shared_mutex _lock;
    };
}

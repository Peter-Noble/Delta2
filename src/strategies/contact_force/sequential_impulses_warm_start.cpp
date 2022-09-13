#include "../../collision_detection/contact.h"
#include "sequential_impulses_warm_start.h"


using namespace Delta2;

SequentialImpulseWarmStart::SequentialImpulseWarmStart() {

}

void SequentialImpulseWarmStart::apply(std::vector<LocalContact>& contacts, std::vector<collision::Contact<double>>& hits) {
    std::shared_lock<std::shared_mutex> guard;

    // for (auto& s : _store) {
    //     printf("(%i, %i): %f\n", s.first.first, s.first.second, s.second);
    // }

    for (int c = 0; c < hits.size(); c++) {
        int lower = std::min(hits[c].p_a->id, hits[c].p_b->id);
        int upper = std::max(hits[c].p_a->id, hits[c].p_b->id);
        std::pair<int, int> key = std::make_pair(lower, upper);

        // printf("(%i, %i)\n", lower, upper);

        if (_store.find(key) != _store.end()) {
            contacts[c].impulse = _store[key] * 0.75;
        }
    }
}

void SequentialImpulseWarmStart::update(std::vector<LocalContact>& contacts, std::vector<collision::Contact<double>>& hits) {
    std::unordered_map<std::pair<int, int>, std::pair<double, int>, index_pair_hash> acc;

    for (int c = 0; c < hits.size(); c++) {
        int lower = std::min(hits[c].p_a->id, hits[c].p_b->id);
        int upper = std::max(hits[c].p_a->id, hits[c].p_b->id);
        std::pair<int, int> key = std::make_pair(lower, upper);

        if (acc.find(key) != acc.end()) {
            acc[key] = std::make_pair(contacts[c].impulse, 1);
        }
        else {
            acc[key].first += contacts[c].impulse;
            acc[key].second++;
        }
    };

    std::unique_lock<std::shared_mutex> guard;

    for (auto& tmp : acc) {
        // printf("Updating: (%i, %i) from %f to %f\n", tmp.first.first, tmp.first.second, _store[tmp.first], tmp.second.first / tmp.second.second);
        _store[tmp.first] = tmp.second.first / tmp.second.second;
    }
}

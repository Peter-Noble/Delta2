#include "contact_detection_hybrid.h"

#include "../../collision_detection/hybrid_dist.h"

using namespace Delta2;
using namespace strategy;

ContactDetectionHybrid::ContactDetectionHybrid() {
    
}

void ContactDetectionHybrid::findContactsBucket(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketHybridDeferred(bucket_pairs, a, b, a.future_state.getTransformation(), b.future_state.getTransformation(), hits, pair_used_out);
}

void ContactDetectionHybrid::findContactsBucketCurrent(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketHybridDeferred(bucket_pairs, a, b, a.current_state.getTransformation(), b.current_state.getTransformation(), hits, pair_used_out);
}

void ContactDetectionHybrid::findContactsBucketLast(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketHybridDeferred(bucket_pairs, a, b, a.last_state.getTransformation(), b.last_state.getTransformation(), hits, pair_used_out);
}

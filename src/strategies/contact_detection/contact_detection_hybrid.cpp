#include "contact_detection_hybrid.h"

#include "../../collision_detection/comparison_dist.h"

using namespace Delta2;
using namespace strategy;

ContactDetectionHybrid::ContactDetectionHybrid() {
    
}

void ContactDetectionHybrid::findContactsBucketComparison(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketComparison<8, double>(bucket_pairs, a, b, hits, pair_used_out);
}

void ContactDetectionHybrid::findContactsBucketComparisonCurrent(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketComparisonCurrent<8, double>(bucket_pairs, a, b, hits, pair_used_out);
}

void ContactDetectionHybrid::findContactsBucketComparisonLast(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) {
    collision::findContactsBucketComparisonLast<8, double>(bucket_pairs, a, b, hits, pair_used_out);
}

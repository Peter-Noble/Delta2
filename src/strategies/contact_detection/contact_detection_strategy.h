#pragma once

#include "../../model/collision_state.h"
#include "../../collision_detection/contact.h"

namespace Delta2 {
    namespace strategy {
        class ContactDetectionStrategy {
        public:
            ContactDetectionStrategy();
            virtual void findContactsBucketComparison(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) = 0;
            virtual void findContactsBucketComparisonCurrent(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) = 0;
        protected:

        };
    }
}

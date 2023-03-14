#pragma once

#include "contact_detection_strategy.h"

namespace Delta2 {
    namespace strategy {
        class ContactDetectionComparison : public ContactDetectionStrategy {
        public:
            ContactDetectionComparison();
            void findContactsBucket(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) override;
            void findContactsBucketCurrent(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) override;
            void findContactsBucketLast(const std::vector<collision::DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<collision::Contact<double>>& hits, std::vector<int>& pair_used_out) override;
        protected:

        };
    }
}

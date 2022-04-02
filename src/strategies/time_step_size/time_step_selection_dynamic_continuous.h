#pragma once

#include "time_step_selection_strategy.h"

#include "../contact_detection_continuous/contact_detection_continuous_strategy.h"

namespace Delta2 {
    namespace strategy {
        class TimeStepSelectionDynamicContinuous : public TimeStepSelectionStrategy {
        public:
            TimeStepSelectionDynamicContinuous(ContactDetectionContinuousStrategy& continuous_contact, common::Options& opt);
            double selectTimeStep(collision::Cluster& cluster) override;
        private:
            double _step_size;
            ContactDetectionContinuousStrategy& _contact_continuous;
        };
    }
}

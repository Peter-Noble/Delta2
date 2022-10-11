#pragma once

#include "PDE_strategy.h"

#include "../contact_detection/contact_detection_strategy.h"
#include "../contact_force/contact_force_strategy.h"
#include "../friction/friction_strategy.h"
#include "../time_step_size/time_step_selection_strategy.h"

namespace Delta2 {
    namespace strategy {
        class PDEExplicit : public PDEStrategy {
        public:
            PDEExplicit(ContactDetectionStrategy& contact_detection, ContactForceStrategy& contact_force, FrictionStrategy& friction, TimeStepSelectionStrategy& time_step, common::Options& opt);
            double selectTimeStep(collision::Cluster& cluster) override;
            bool step(collision::Cluster& cluster, bool allow_fail=true) override;
            void printType() override {
                printf("Explicit\n");
            };
        protected:
            ContactDetectionStrategy& _contact_detection;
            ContactForceStrategy& _contact_force;
            FrictionStrategy& _friction;
            TimeStepSelectionStrategy& _time_step;
        };
    }
}

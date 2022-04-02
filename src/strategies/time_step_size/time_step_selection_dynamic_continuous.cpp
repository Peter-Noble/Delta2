#include "time_step_selection_dynamic_continuous.h"

#include "../../collision_detection/separate_clusters.h"

using namespace Delta2;
using namespace strategy;


TimeStepSelectionDynamicContinuous::TimeStepSelectionDynamicContinuous(ContactDetectionContinuousStrategy& contact_continuous,
                                                                       common::Options& opt) :
                                                                       TimeStepSelectionStrategy(opt),
                                                                       _step_size(opt.time_step_size),
                                                                       _contact_continuous(contact_continuous) {

}

double TimeStepSelectionDynamicContinuous::selectTimeStep(collision::Cluster& cluster) {
    collision::fineCollisionClustersWithTimeStepSelection(cluster);
    cluster.step_size = std::min(cluster.step_size, (double)_opt.time_step_size);
    return cluster.step_size;
}

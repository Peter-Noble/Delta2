#include "time_step_selection_strategy.h"

using namespace Delta2;
using namespace strategy;

TimeStepSelectionStrategy::TimeStepSelectionStrategy(common::Options& opt) :
                                                     _opt(opt) {

}

void TimeStepSelectionStrategy::init(model::ParticleHandler& ps) {
    for (Particle* p : ps) {
        p->last_time_step_size = _opt.time_step_size;
        p->last_state = p->current_state;
        p->future_state = p->current_state;
    }
}

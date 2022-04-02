#include "broad_phase_strategy.h"

using namespace Delta2::strategy;


BroadPhaseStrategy::BroadPhaseStrategy(PDEStrategy& local_pde,
                                       common::Options& opt) :
                                       _local_pde(local_pde),
                                       _opt(opt)
{

}

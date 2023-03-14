#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../globals.h"
#include "../common/cli_options.h"
#include <mutex>
#include <Eigen/Dense>

using namespace Delta2;

std::mutex globals::contact_draws_lock;
std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> globals::contact_draws;
Delta2::common::Options globals::opt;
Delta2::common::IttHandles globals::itt_handles;
Delta2::common::Logger globals::logger;

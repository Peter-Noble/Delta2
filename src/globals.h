#pragma once

#include <vector>
#include <Eigen/Dense>
#include <mutex>
#include "./common/cli_options.h"
#include "./common/itt.h"
#include "./common/logger.h"

namespace Delta2 {
    namespace globals {
        struct Step {
            int depth;
            float size;
            void print() {
                printf("%i, %f\n", depth, size);
            }
        };

        struct Rollback {
            int depth; // -1 if rendezvous
            float parent_step;
            float rollback_by;
            void print() {
                printf("%i, %f, %f\n", depth, parent_step, rollback_by);
            }
        };

        struct ClusterInfo {
            int num_particles;
            int num_interactions;
            int num_contacts;
            void print() {
                printf("%i, %i, %i\n", num_particles, num_interactions, num_contacts);
            }
        };

        struct Stats {
            std::vector<Step> time_steps;
            std::vector<Rollback> rollbacks;
            std::vector<ClusterInfo> clusters;
            std::mutex m;

            void print() {
                printf("==== STATS ====\n");
                printf("==== time steps ====\n");
                for (Step& s : time_steps) {
                    s.print();
                }
                printf("==== rollbacks ====\n");
                for (Rollback& r : rollbacks) {
                    r.print();
                }
                printf("==== clusters ====\n");
                for (ClusterInfo& c : clusters) {
                    c.print();
                }
            }
            void addTimestep(Step s) {
                std::lock_guard lock(m);
                time_steps.push_back(s);
            }
            void addRollback(Rollback s) {
                std::lock_guard lock(m);
                rollbacks.push_back(s);
            }
            void addCluster(ClusterInfo s) {
                std::lock_guard lock(m);
                clusters.push_back(s);
            }
        };
    }
}

namespace Delta2 {
    namespace globals {
        extern std::mutex contact_draws_lock;
        extern std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> contact_draws;
        extern Delta2::common::Options opt;
        extern Delta2::common::IttHandles itt_handles;
        extern Delta2::common::Logger logger;
        extern Delta2::globals::Stats stats;
    }
}

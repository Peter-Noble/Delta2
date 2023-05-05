#pragma once

#include "colour_hits.h"
#include "broad_phase.h"
#include "../model/particle.h"
#include "../model/particle_handler.h"
#include "contact.h"
#include "contact_state.h"
#include "cluster.h"

#include "tbb/tbb.h"

namespace Delta2 {
    namespace collision {
        std::vector<uint32_t> simpleParallelConnectedComponents(collision::BroadPhaseCollisions& broad_phase, model::ParticleHandler& particles) {
            int num_particles = particles.size();
            std::vector<uint32_t> labels;
            labels.resize(num_particles);

            std::vector<Contact<double>> hits;
            for (BroadPhaseCollision& b : broad_phase) {
                int a_id = particles.getLocalID(b.A.first); // local geo id
                int b_id = particles.getLocalID(b.B.first);

                if (!particles[a_id].is_static && !particles[b_id].is_static) {
                    hits.emplace_back(
                        Eigen::Vector3d({0, 0, 0}),
                        Eigen::Vector3d({0, 0, 0}),
                        0,
                        0,
                        0,
                        0,
                        particles[a_id],
                        particles[b_id]
                    );
                }
            }

            std::vector<std::vector<ContactBundle>> colours = collision::colour_hits(particles, hits);

            // For each particle set current label
            tbb::parallel_for(tbb::blocked_range<int>(0, num_particles), [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i < r.end(); i++) {
                    labels[i] = i;
                }
            });

            // For each particle look at neightbours to select minimum id
            // const int parallel_threshold = globals::opt.sequential_parallel_threshold;
            const int parallel_individual_colour_threshold = 500;
            const int parallel_grain_size = 100;

            bool no_change = false;

            while (no_change) {
                no_change = true;
                for (int colour = 0; colour < colours.size(); colour++) {
                    if (colours[colour].size() > parallel_individual_colour_threshold) {
                        int max_concurrency = std::min(tbb::info::default_concurrency(), (int)colours[colour].size() / parallel_grain_size);
                        tbb::task_arena arena(max_concurrency);
                        arena.execute([&] {
                            tbb::parallel_for(tbb::blocked_range<int>(0, colours[colour].size(), parallel_grain_size),
                                            [&](tbb::blocked_range<int> r) {
                                for (int i = r.begin(); i < r.end(); i++) {
                                    uint32_t label_a = labels[colours[colour][i].lower];
                                    uint32_t label_b = labels[colours[colour][i].upper];
                                    if (label_a != label_b) {
                                        no_change = false;
                                        if (label_a < label_b) {
                                            labels[colours[colour][i].upper] = label_a;
                                        }
                                        else {
                                            labels[colours[colour][i].lower] = label_b;
                                        }
                                    }
                                }
                            });
                        });
                    }
                    else {
                        for (int b = 0; b < colours[colour].size(); b++) {
                            uint32_t label_a = labels[colours[colour][b].lower];
                            uint32_t label_b = labels[colours[colour][b].upper];
                            if (label_a != label_b) {
                                no_change = false;
                                if (label_a < label_b) {
                                    labels[colours[colour][b].upper] = label_a;
                                }
                                else {
                                    labels[colours[colour][b].lower] = label_b;
                                }
                            }
                        }
                    }
                }
            }
            return labels;
        }
    }
}

#include "parallel_connected_components.h"
#include "../globals.h"

using namespace Delta2;

std::vector<uint32_t> collision::simpleParallelConnectedComponents(collision::BroadPhaseCollisions& broad_phase, model::ParticleHandler& particles) {
    __itt_string_handle* init_hits_task = __itt_string_handle_create("parallel connected init hits");

    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, init_hits_task);
    int num_particles = particles.size();
    std::vector<uint32_t> labels;
    labels.resize(num_particles);

    // std::vector<Contact<double>> hits;
    // for (BroadPhaseCollision& b : broad_phase) {
    //     int a_id = particles.getLocalID(b.A.first); // local geo id
    //     int b_id = particles.getLocalID(b.B.first);

    //     if (!particles[a_id].is_static && !particles[b_id].is_static) {
    //         hits.emplace_back(
    //             Eigen::Vector3d({0, 0, 0}),
    //             Eigen::Vector3d({0, 0, 0}),
    //             0,
    //             0,
    //             0,
    //             0,
    //             particles[a_id],
    //             particles[b_id]
    //         );
    //     }
    // }
    std::vector<ContactBundle> bundles;
    for (BroadPhaseCollision& b : broad_phase) {
        const int a_id = particles.getLocalID(b.A.first); // local geo id
        const int b_id = particles.getLocalID(b.B.first);
        // const int lower = std::min(a_id, b_id);
        // const int upper = std::max(a_id, b_id);
        const int lower = a_id;
        const int upper = b_id;

        if (!particles[a_id].is_static && !particles[b_id].is_static) {
            bundles.emplace_back();
            const int b = bundles.size() - 1;
            bundles[b].lower = lower;
            bundles[b].upper = upper;
            bundles[b].normal_average = {0, 0, 0};
            bundles[b].tangent_average = {0, 0, 0};
        }
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    globals::logger.printf(1, "Colour hits\n");
    // std::vector<std::vector<ContactBundle>> colours = collision::colour_hits(particles, hits);
    std::vector<std::vector<ContactBundle>> colours = collision::colour_hits(particles, bundles);
    globals::logger.printf(1, "Colour hits done\n");

    // For each particle set current label
    for (int i = 0; i < num_particles; i++) {
        labels[i] = i;
    }
    // tbb::parallel_for(tbb::blocked_range<int>(0, num_particles), [&](tbb::blocked_range<int> r) {
    //     for (int i = r.begin(); i < r.end(); i++) {
    //         labels[i] = i;
    //     }
    // });

    // For each particle look at neightbours to select minimum id
    // const int parallel_threshold = globals::opt.sequential_parallel_threshold;
    const int parallel_individual_colour_threshold = 200;
    const int parallel_grain_size = 100;

    bool no_change = false;

    globals::logger.printf(1, "Loop begin\n");
    while (!no_change) {
        globals::logger.printf(1, "Inner loop\n");
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
    globals::logger.printf(1, "Loop end\n");
    return labels;
}
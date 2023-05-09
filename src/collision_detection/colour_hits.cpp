#include "colour_hits.h"
#include "../globals.h"

using namespace Delta2;
using namespace collision;

// A hash function used to hash a pair of any kind
struct HashPair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

std::vector<ContactBundle> Delta2::collision::bundle_hits(model::ParticleHandler& particles, std::vector<Contact<double>>& hits) {
    __itt_string_handle* bundle_task = __itt_string_handle_create("bundle_hits");
    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, bundle_task);

    std::vector<ContactBundle> bundles;

    std::unordered_map<std::pair<int, int>, int, HashPair> bundle_lookup;

    int last_lower = -1;
    int last_upper = -1;
    
    // std::sort(hits.begin(), hits.end(), [](const Contact<double>& a, const Contact<double>& b) {
    //     return a.p_a->id < b.p_a->id || (a.p_a->id == b.p_a->id && a.p_b->id < b.p_b->id);
    // });  // Even for very small numbers this is extraordinarily slow for some reason.  Deep copying data?
    // Instead do sort while combining into bundles

    int b = -1;
    for (int c = 0; c < hits.size(); c++) {
        const int lower = std::min(hits[c].p_a->id, hits[c].p_b->id);
        const int upper = std::max(hits[c].p_a->id, hits[c].p_b->id);

        if (lower != last_lower || upper != last_upper) {
            bool found = false;
            const int local_lower = particles.getLocalID(lower);
            const int local_upper = particles.getLocalID(upper);
            // auto key = std::make_pair(local_lower, local_upper);
            // auto search = bundle_lookup.find(key);
            // if (search == bundle_lookup.end()) {
            //     bundles.emplace_back();
            //     b = bundles.size() - 1;
            //     bundles[b].lower = local_lower;
            //     bundles[b].upper = local_upper;
            //     bundles[b].normal_average = {0, 0, 0};
            //     bundles[b].tangent_average = {0, 0, 0};
                
            //     bundle_lookup.insert({key, b});
            // } else {
            //     b = search->second;
            // }
            for (int b_old = 0; b_old < bundles.size(); b_old++) {
                if (bundles[b_old].lower == local_lower && bundles[b_old].upper == local_upper) {
                    found = true;
                    b = b_old;
                    break;
                }
            }
            if (!found) {
                bundles.emplace_back();
                b = bundles.size() - 1;
                bundles[b].lower = local_lower;
                bundles[b].upper = local_upper;
                bundles[b].normal_average = {0, 0, 0};
                bundles[b].tangent_average = {0, 0, 0};
            }

            last_lower = lower;
            last_upper = upper;
        }

        bundles[b].hits.push_back(c);
        assert(c < hits.size());
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    return bundles;
}
std::vector<std::vector<ContactBundle>> Delta2::collision::colour_hits(model::ParticleHandler& particles, std::vector<Contact<double>>& hits) {
    std::vector<ContactBundle> bundles = Delta2::collision::bundle_hits(particles, hits);
    return colour_hits(particles, bundles);
}


std::vector<std::vector<ContactBundle>> Delta2::collision::colour_hits(model::ParticleHandler& particles, std::vector<ContactBundle>& bundles) {
    __itt_string_handle* adjacency_task = __itt_string_handle_create("Colour hits - adjacency");
    __itt_string_handle* max_degree_task = __itt_string_handle_create("Colour hits - max_degree");
    __itt_string_handle* init_avail_task = __itt_string_handle_create("Colour hits - init_avail");
    __itt_string_handle* main_loop_task = __itt_string_handle_create("Colour hits - main_loop");

    // Returns a vector for each colour used containing the indices of the hits with that colour groups by batch (same pair of involved particles)

    std::vector<int> particle_adj[particles.size()];
    
    int max_degree = 0;

    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, adjacency_task);
    int colour[bundles.size()];
    for (int b = 0; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        if (!particles[a_id].is_static) {
            particle_adj[a_id].push_back(b);
        }
        if (!particles[b_id].is_static) {
            particle_adj[b_id].push_back(b);
        }

        colour[b] = -1;
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, max_degree_task);
    for (int b = 0; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        const int degree = particle_adj[a_id].size() + particle_adj[b_id].size() - 2;

        max_degree = std::max(max_degree, degree);
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    int max_used = 0;

    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, init_avail_task);
    std::vector<std::vector<ContactBundle>> result;
    bool avail[max_degree+1]; // Greedy algorithm guarentees degree+1 maximum colours
    for (int i = 0; i < max_degree+1; i++) {
        avail[i] = true;
        result.push_back({});
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    colour[0] = 0;
    result[0].push_back(bundles[0]);

    __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, main_loop_task);
    for (int b = 1; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        if (!particles[a_id].is_static) {
            for (int a : particle_adj[a_id]) {
                if (colour[a] != -1) {
                    avail[colour[a]] = false;
                }
            }
        }
        if (!particles[b_id].is_static) {
            for (int a : particle_adj[b_id]) {
                if (colour[a] != -1) {
                    avail[colour[a]] = false;
                }
            }
        }

        int found = 0;
        while (!avail[found]) found++;

        colour[b] = found;
        result[found].push_back(bundles[b]);
        max_used = std::max(max_used, found + 1);

        for (int a : particle_adj[a_id]) {
            avail[colour[a]] = true;
        }
        for (int a : particle_adj[b_id]) {
            avail[colour[a]] = true;
        }
    }
    __itt_task_end(globals::itt_handles.detailed_domain);

    result.resize(max_used);

    globals::logger.printf(3, "Colouring graph of size %i bundles with %i colours from max degree %i\n", bundles.size(), result.size(), max_degree);

    return result;
}

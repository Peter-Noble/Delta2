#include "colour_hits.h"
#include "../globals.h"

using namespace Delta2;
using namespace collision;

std::vector<ContactBundle> Delta2::collision::bundle_hits(Cluster& cluster, std::vector<Contact<double>>& hits) {
    std::vector<ContactBundle> bundles;

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
            const int local_lower = cluster.particles.getLocalID(lower);
            const int local_upper = cluster.particles.getLocalID(upper);
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

    return bundles;
}

std::vector<std::vector<ContactBundle>> Delta2::collision::colour_hits(Cluster& cluster, std::vector<Contact<double>>& hits) {
    // Returns a vector for each colour used containing the indices of the hits with that colour groups by batch (same pair of involved particles)
    std::vector<ContactBundle> bundles = Delta2::collision::bundle_hits(cluster, hits);

    std::vector<int> particle_adj[cluster.particles.size()];
    
    int max_degree = 0;

    int colour[bundles.size()];
    for (int b = 0; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        if (!cluster.particles[a_id].is_static) {
            particle_adj[a_id].push_back(b);
        }
        if (!cluster.particles[b_id].is_static) {
            particle_adj[b_id].push_back(b);
        }

        colour[b] = -1;
    }

    for (int b = 0; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        const int degree = particle_adj[a_id].size() + particle_adj[b_id].size() - 2;

        max_degree = std::max(max_degree, degree);
    }

    int max_used = 0;

    std::vector<std::vector<ContactBundle>> result;
    bool avail[max_degree+1]; // Greedy algorithm guarentees degree+1 maximum colours
    for (int i = 0; i < max_degree+1; i++) {
        avail[i] = true;
        result.push_back({});
    }

    colour[0] = 0;
    result[0].push_back(bundles[0]);

    for (int b = 1; b < bundles.size(); b++) {
        const int a_id = bundles[b].lower;
        const int b_id = bundles[b].upper;

        if (!cluster.particles[a_id].is_static) {
            for (int a : particle_adj[a_id]) {
                if (colour[a] != -1) {
                    avail[colour[a]] = false;
                }
            }
        }
        if (!cluster.particles[b_id].is_static) {
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

    result.resize(max_used);

    globals::logger.printf(3, "Colouring graph of size %i (%i bundles) with %i colours from max degree %i\n", hits.size(), bundles.size(), result.size(), max_degree);

    return result;
}

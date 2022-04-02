#include "../model/surrogate_components.h"
#include "../model/collision_state.h"
#include "../collision_detection/comparison_dist.h"

using namespace Delta2::model;

std::vector<Delta2::collision::Contact<double>> compare_surrogate_tree(Delta2::Particle a, Delta2::Particle b, Delta2::collision::ActiveSet active_a, Delta2::collision::ActiveSet active_b, std::vector<Delta2::collision::Contact<double>>& intermediate_hits, bool& more_layers_to_expand, bool& intermediate_hit) {
    std::vector<Delta2::collision::DeferredCompare> dcs;
    std::vector<Delta2::collision::Contact<double>> hits;
    std::vector<Delta2::collision::Contact<double>> result;
    std::vector<int> used;
    SurrogateTree surrogate_a = a.mesh->getSurrogateTree();
    SurrogateTree surrogate_b = b.mesh->getSurrogateTree();
    for (const Delta2::collision::NodeState& state_a : active_a.states) {
        Node node_a = surrogate_a.getNode(state_a.node_id);
        for (const Delta2::collision::NodeState& state_b : active_b.states) {
            Node node_b = surrogate_b.getNode(state_b.node_id);
            dcs.emplace_back(node_a, node_b);
        }
    }
    Delta2::collision::findContactsBucketComparison<8, double>(dcs, a, b, hits, used);
    intermediate_hit = false;
    int last_used = 0;
    for (int dc_i = 0; dc_i < dcs.size(); dc_i++) {
        const Delta2::collision::DeferredCompare& dc = dcs[dc_i];
        int state_a_id = dc_i / active_b.states.size();
        int state_b_id = dc_i % active_b.states.size();
        if (used[dc_i] > last_used) {
            active_a.states[state_a_id].expand = true;
            active_b.states[state_b_id].expand = true;

            if (dcs[dc_i].a.is_inner && dcs[dc_i].b.is_inner) {
                result.insert(result.end(), hits.begin() + last_used, hits.begin() + used[dc_i]);
            } 
            else if (active_a.states[state_a_id].newly_expanded || active_b.states[state_b_id].newly_expanded) {
                intermediate_hits.insert(intermediate_hits.end(), hits.begin() + last_used, hits.begin() + used[dc_i]);
                intermediate_hit = true;
            }
        } else {
            active_a.states[state_a_id].reduce = true;
            active_b.states[state_b_id].reduce = true;
        }
        last_used = used[dc_i];
    }
    return hits;
}

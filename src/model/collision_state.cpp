#include "collision_state.h"
#include "surrogate_components.h"

using namespace Delta2::collision;

PairState::PairState(Particle* A, Particle* B) {
    /*a_current_rotation = A->rotation;
	a_current_angular_momentum = A->angular_momentum;
	a_current_angular_velocity = A->angular_velocity;
	a_current_translation = A->translation;
	a_current_velocity = A->linear_velocity;

	b_current_rotation = B->rotation;
	b_current_angular_momentum = B->angular_momentum;
	b_current_angular_velocity = B->angular_velocity;
	b_current_translation = B->translation;
	b_current_velocity = B->linear_velocity;*/

	a_last_force = { 0, 0, 0 };
	a_avg_force = { 0, 0, 0 };
	a_last_update_force = { 0, 0, 0 };
	a_last_torque = { 0, 0, 0 };
	a_avg_torque = { 0, 0, 0 };
	a_last_update_torque = { 0, 0, 0 };
	a_force = { 0, 0, 0 };
	a_torque = { 0, 0, 0 };

	b_last_force = { 0, 0, 0 };
	b_avg_force = { 0, 0, 0 };
	b_last_update_force = { 0, 0, 0 };
	b_last_torque = { 0, 0, 0 };
	b_avg_torque = { 0, 0, 0 };
	b_last_update_torque = { 0, 0, 0 };
	b_force = { 0, 0, 0 };
	b_torque = { 0, 0, 0 };

	/*a_a = { 0, 0, 0 };
	a_t = { 0, 0, 0 };
	b_a = { 0, 0, 0 };
	b_t = { 0, 0, 0 };*/
}

ActiveSet updateActiveSet(const ActiveSet& S, Delta2::model::SurrogateTree* T) {
	ActiveSet result;
	std::unordered_map<int, std::vector<int>> parents_to_reduce;
	for (const NodeState& s : S.states) {
		const Delta2::model::Node& node = T->getNode(s.node_id);
		if (s.expand) {
			for (int c = 0; c < node.num_children; c++) {
				NodeState n;
				n.node_id = c + node.child_id_first;
				n.expand = false;
				n.reduce = false;
				n.newly_expanded = true;
				result.states.push_back(n);
			}
		}
		else if (s.reduce) {
			parents_to_reduce[node.parent_id].push_back(s.node_id);
		}
	}
	for (const auto& [node_id, children] : parents_to_reduce) {
		const Delta2::model::Node& parent = T->getNode(node_id);
		if (children.size() == parent.num_children) {
			NodeState n;
			n.node_id = node_id;
			n.expand = false;
			n.reduce = false;
			n.newly_expanded = false;
			result.states.push_back(n);
		}
		else {
			for (int c : children) {
				NodeState n;
				n.node_id = c;
				n.expand = false;
				n.reduce = false;
				n.newly_expanded = false;
				result.states.push_back(n);
			}
		}
	}
	return result;
}

void updateActiveSet(const ActiveSet& S) {
	
}

DeferredCompare::DeferredCompare(const Delta2::model::Node& first, const Delta2::model::Node& second) : a(first), b(second) {};

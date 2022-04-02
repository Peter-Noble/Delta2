#include "contact_state.h"
#include "../model/collision_state.h"

using namespace Delta2;

collision::SeparationWitness::SeparationWitness(uint32_t a_n, uint32_t b_n, const Eigen::Vector3d& orig, const Eigen::Vector3d& norm, bool a_def) {
    a_node = a_n;
    b_node = b_n;
    origin = orig;
    normal = norm;
    a_defining = a_def;
}

bool collision::SeparationWitness::check(Particle& a, Particle& b) {
    uint32_t d_node = a_defining ? a_node : b_node;
    uint32_t o_node = a_defining ? b_node : a_node;

    Particle& d = a_defining ? a : b;
    Particle& o = a_defining ? b : a;
    
    const model::Node& d_n = d.mesh->getSurrogateTree().getNode(d_node);
    const model::Node& o_n = o.mesh->getSurrogateTree().getNode(d_node);

    std::shared_ptr<model::Bucket> d_bucket = d.mesh->getSurrogateTree().getBucket(d_n.bucket_id);
    std::shared_ptr<model::Bucket> o_bucket = o.mesh->getSurrogateTree().getBucket(o_n.bucket_id);

    // std::vector<Eigen::Vector3d> defining_verts = d_bucket->getVertices(); // We don't have to test these.  By definition they're on the right side of the plane.
    std::vector<Eigen::Vector3d> other_verts = o_bucket->getVertices();

    Eigen::Vector3d orig_current = common::transform(origin, d.current_state.getTransformation());
    Eigen::Vector3d normal_current = common::transform(normal, d.current_state.getRotation().toRotationMatrix());
    normal_current.normalize();

    Eigen::Vector3d orig_future = common::transform(origin, d.future_state.getTransformation());
    Eigen::Vector3d normal_future = common::transform(normal, d.future_state.getRotation().toRotationMatrix());
    normal_future.normalize();

    bool separated = true;

    for (const Eigen::Vector3d& v : other_verts) {
        if ((common::transform(v, o.current_state.getTransformation()) - orig_current).dot(normal_current) < eps) {
            separated = false;
            break;
        }
        if ((common::transform(v, o.future_state.getTransformation()) - orig_future).dot(normal_future) < eps) {
            separated = false;
            break;
        }
    }
    // TODO we might gain some valuable information here about which vertices might produce contacts.

    return separated;    
}

collision::DeferredCompare collision::SeparationWitness::toDeferredCompare(Particle& a, Particle& b) {
    collision::DeferredCompare dc(a.mesh->getSurrogateTree().getNode(a_node), b.mesh->getSurrogateTree().getNode(b_node));
    return dc;
}


collision::ContactState::ContactState() {

}

collision::ContactState::ContactState(collision::DeferredCompare dc) {
    noWitnesses.push_back(dc);
}

bool collision::ContactState::attemptWitness(Particle& a, Particle& b, const DeferredCompare& dc) {
    std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
    std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

    double eps = a_bucket->getMaxEps() + b_bucket->getMaxEps(); // TODO use the individual triangle eps (this may be a combination of different triangle epss)

    Eigen::Vector3d plane_normal, plane_origin;

    bool direction_chosen = false;
    for (const common::Triangle<double>& t : a_bucket->getConvexHull()) {
        auto [plane_normal_current, plane_origin_current] = t.transformed(a.current_state.getTransformation()).toNormalRep();
        auto [plane_normal_future, plane_origin_future] = t.transformed(a.future_state.getTransformation()).toNormalRep();

        bool valid = true;
        
        for (const Eigen::Vector3d v : b_bucket->getVertices()) {
            double dist = (common::transform(v, b.current_state.getTransformation()) - plane_origin_current).dot(plane_normal_current);
            if (!direction_chosen) {
                direction_chosen = true;
                if (dist < 0.0) {
                    plane_normal_current *= -1;
                    plane_normal_future *= -1;
                    dist *= -1;
                }
            }
            if (dist < eps) {
                valid = false;
                break;
            }
            dist = (common::transform(v, b.future_state.getTransformation()) - plane_origin_future).dot(plane_normal_future);

        }

        if (valid) {
            auto plane_local = t.toNormalRep();
            plane_normal = plane_local.first;
            plane_origin = plane_local.second;
            collision::SeparationWitness sw(dc.a.node_id, dc.b.node_id, plane_origin, plane_normal, true);
            witnesses.push_back(sw);
            return true;
        }
    }

    direction_chosen = false;
    for (const common::Triangle<double>& t : b_bucket->getConvexHull()) {
        auto [plane_normal_current, plane_origin_current] = t.transformed(b.current_state.getTransformation()).toNormalRep();
        auto [plane_normal_future, plane_origin_future] = t.transformed(b.future_state.getTransformation()).toNormalRep();

        bool valid = true;
        
        for (const Eigen::Vector3d v : a_bucket->getVertices()) {
            double dist = (common::transform(v, a.current_state.getTransformation()) - plane_origin_current).dot(plane_normal_current);
            if (!direction_chosen) {
                direction_chosen = true;
                if (dist < 0.0) {
                    plane_normal_current *= -1;
                    plane_normal_future *= -1;
                    dist *= -1;
                }
            }
            if (dist < eps) {
                valid = false;
                break;
            }
            dist = (common::transform(v, a.future_state.getTransformation()) - plane_origin_future).dot(plane_normal_future);

        }

        if (valid) {
            auto plane_local = t.toNormalRep();
            plane_normal = plane_local.first;
            plane_origin = plane_local.second;
            collision::SeparationWitness sw(dc.a.node_id, dc.b.node_id, plane_origin, plane_normal, false);
            witnesses.push_back(sw);
            return true;
        }
    }

    // No witness was found
    noWitnesses.push_back(dc);
    return false;
}

collision::ContactStateCache::ContactStateCache() {
    
}

collision::ContactState& collision::ContactStateCache::getContactState(Particle& a, Particle& b) {
    std::lock_guard<std::mutex> guard(_lock);
    auto key = std::make_pair(a.id, b.id);
    auto search = states.find(key);
    if (search == states.end()) {
        collision::DeferredCompare dc(a.mesh->getSurrogateTree().getNode(0), b.mesh->getSurrogateTree().getNode(0));
        collision::ContactState cs(dc);
        auto iter = states.insert({key, cs});
        return iter.first->second;
    } else {
        return search->second;
    }
}

#include <vector>

namespace Delta2 {
    namespace collision {
        using SurfacePoint = std::pair<unsigned,unsigned>;  // geom id, prim id
        using BroadPhaseCollision = std::pair<SurfacePoint, SurfacePoint>;
        using BroadPhaseCollisions = std::vector<BroadPhaseCollision>;
    }
}

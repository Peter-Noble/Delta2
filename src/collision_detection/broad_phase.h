#pragma once

#include <vector>

namespace Delta2 {
    namespace collision {
        using SurfacePoint = std::pair<unsigned,unsigned>;  // geom id, prim id
        // using BroadPhaseCollision = std::pair<SurfacePoint, SurfacePoint>;
        struct BroadPhaseCollision {
            SurfacePoint A;
            SurfacePoint B;
            float min_toc;
            float target_toc;
        };
        using BroadPhaseCollisions = std::vector<BroadPhaseCollision>;
    }
}

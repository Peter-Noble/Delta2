#include "../testing/unit_test.h"
#include "particle.h"
#include "../common/primitive_geo.h"
#include "../testing/approx_bbox.h"

TEST_CASE("Particle creation", "[grp_model][particle]") {
    Delta2::common::Options opt;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Delta2::common::cube(1, V, F);

    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    Delta2::Particle A(M, 1.0, 1.0, 0.01);

    REQUIRE( A.getMass() == Approx(8.0) );
    REQUIRE_THAT( A.mesh->getAABB(), BBoxEqual(embree::BBox3fa({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0})) );
    float radius = sqrt(3);
    REQUIRE_THAT( A.mesh->getSphereBounds(), BBoxEqual(embree::BBox3fa({-radius, -radius, -radius}, {radius, radius, radius})) );
}

TEST_CASE("Particle point velocity", "[grp_model][particle]") {
    Delta2::common::Options opt;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Delta2::common::cube(1, V, F);

    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    Delta2::Particle A(M, 1.0, 1.0, 0.01);

    A.future_state.setAngular({1.0, 0.0, 0.0});
    Eigen::Vector3d v = A.futurePointVelocity({0.0, 0.0, -1.0});
    REQUIRE( v.normalized().dot(Eigen::Vector3d({0.0, 1.0, 0.0})) == Approx(1.0) );

    A.future_state.setVelocity({0.0, 1.0, 0.0});
    v = A.futurePointVelocity({0.0, 0.0, -1.0});
    REQUIRE( v.normalized().dot(Eigen::Vector3d({0.0, 1.0, 0.0})) == Approx(1.0) );

}

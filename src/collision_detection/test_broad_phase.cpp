#include "../testing/unit_test.h"
#include "broad_phase_embree.h"
#include "../common/primitive_geo.h"

TEST_CASE("Broad phase Embree stationary", "[grp_collision_detection][broad_phase]") {
    std::vector<Delta2::Particle> particles;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Delta2::common::cube(V, F);

    Delta2::common::Options opt;

    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    for (int i = 0; i < 10; i++) {
        Delta2::Particle P(M, 1.0, 1.0, 0.01);
        P.current_state.setTranslation({2.0*i, 0.0, 0.0});
        particles.push_back(P);
    }

    Delta2::collision::BroadPhaseCollisions C = Delta2::collision::broadPhaseEmbree(particles, 1.0);

    REQUIRE( C.size() == 9 );
}


TEST_CASE("Broad phase Embree moving", "[grp_collision_detection][broad_phase]") {
    std::vector<Delta2::Particle> particles;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Delta2::common::cube(V, F);

    Delta2::common::Options opt;

    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, opt));
    
    {
        Delta2::Particle P(M, 1.0, 1.0, 0.01);
        P.current_state.setTranslation({0.0, 0.0, 0.0});
        P.current_state.setVelocity({1.0, 0.0, 0.0});
        particles.push_back(P);
    }
    {
        Delta2::Particle P(M, 1.0, 1.0, 0.01);
        P.current_state.setTranslation({4.0, 0.0, 0.0});
        P.current_state.setVelocity({-1.0, 0.0, 0.0});
        particles.push_back(P);
    }

    Delta2::collision::BroadPhaseCollisions C = Delta2::collision::broadPhaseEmbree(particles, 0.1);

    REQUIRE( C.size() == 0 );
    
    C = Delta2::collision::broadPhaseEmbree(particles, 0.5);

    REQUIRE( C.size() == 1 );
}

TEST_CASE("Broad phase Embree static", "[grp_collision_detection][broad_phase]") {
    std::vector<Delta2::Particle> particles;
    Delta2::common::Options opt;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Delta2::common::cube(V, F);
    std::shared_ptr<Delta2::MeshData> m_cube(new Delta2::MeshData(V, F, opt));

    Delta2::common::plane(10, V, F);
    std::shared_ptr<Delta2::MeshData> m_plane(new Delta2::MeshData(V, F, opt, true));
    
    {
        Delta2::Particle P(m_plane, 1.0, 1.0, 0.01);
        P.is_static = true;
        particles.push_back(P);
    }
    {
        Delta2::Particle P(m_cube, 1.0, 1.0, 0.01);
        P.current_state.setTranslation({0.0, 0.0, 2.0});
        P.current_state.setVelocity({0.0, 0.0, -1.0});
        particles.push_back(P);
    }

    Delta2::collision::BroadPhaseCollisions C = Delta2::collision::broadPhaseEmbree(particles, 0.1);

    REQUIRE( C.size() == 0 );
    
    C = Delta2::collision::broadPhaseEmbree(particles, 1.0);

    REQUIRE( C.size() == 1 );

    {
        Delta2::Particle P(m_plane, 1.0, 1.0, 0.01);
        P.is_static = true;
        particles.push_back(P);
    }

    C = Delta2::collision::broadPhaseEmbree(particles, 1.0);

    REQUIRE( C.size() == 2 );
}

TEST_CASE("Broad phase embree empty", "[grp_collision_detection][broad_phase]") {
    std::vector<Delta2::Particle> particles;

    Delta2::collision::BroadPhaseCollisions C = Delta2::collision::broadPhaseEmbree(particles, 1.0);

    REQUIRE( C.size() == 0 );
}
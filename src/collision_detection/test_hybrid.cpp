#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "penalty_dist.h"
#include "hybrid_dist.h"
#include "../common/primitive_geo.h"
#include "../collision_detection/contact.h"
#include "../globals.h"

TEST_CASE("Triangle-triangle hybrid vert-face", "[grp_collision_detection][hybrid]") {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    // Eigen::Vector3f a_rot({0.5235, 0.3491, 0.0}); // Has a problem with rotation order
    Eigen::Quaternionf a_rot({0.951261, 0.254839, 0.1677497, -0.0449395});
    Eigen::Vector3f a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);
    std::vector<Delta2::common::Triangle<float>> As = {A};

    // Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternionf b_rot({0.9952, -0.087056, 0.0038, 0.04347});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};

    std::vector<float> max_error = { 0.01 };

    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_hybridDists(As, Bs, max_error);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector3f({-0.342, -0.732, 0.326})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector3f({-0.339, -0.767, 0.129})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector3f({-0.342, -0.732, 0.326}) - Eigen::Vector3f({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
}

TEST_CASE("Triangle-triangle hybrid edge-edge", "[grp_collision_detection][hybrid]") {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    // Eigen::Vector3f a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Quaternionf a_rot({0.744, 0.498, 0.368, 0.252});
    Eigen::Vector3f a_pos({-1.0284, -1.0942, 1.1695});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);
    std::vector<Delta2::common::Triangle<float>> As = {A};

    // Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternionf b_rot({0.9952, -0.087056, -0.0038, 0.04347});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};

    std::vector<float> max_error = { 0.01 };
    std::vector<bool> failed = { false };

    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_hybridDists(As, Bs, max_error);

    REQUIRE( result.size() == 1 );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector3f({-0.6545, -0.9378, 0.2656}) - Eigen::Vector3f({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector3f({-0.6545, -0.9378, 0.2656})).margin(0.05) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector3f({-0.5089, -0.8695, 0.1449})).margin(0.05) );
}

TEST_CASE("Triangle-triangle hybrid intersecting", "[grp_collision_detection][hybrid]") {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    // Eigen::Vector3f a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Quaternionf a_rot({0.744, 0.498, 0.368, 0.252});
    Eigen::Vector3f a_pos({0.052572, -0.153374, 1.0443});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);
    std::vector<Delta2::common::Triangle<float>> As = {A};

    // Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternionf b_rot({0.9952, -0.087056, -0.0038, 0.04347});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};

    std::vector<float> max_error = { 0.01 };
    std::vector<bool> failed = { false };

    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_hybridDists(As, Bs, max_error);

    REQUIRE( result.size() == 1 );
    REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
}

TEST_CASE("Triangle-triangle hybrid bucket", "[grp_collision_detection][hybrid]") {
    Eigen::Vector3d a({-1, 0, 0});
    Eigen::Vector3d b({1, 0, 0});
    Eigen::Vector3d c({0, -2, 0});
    Delta2::common::Triangle<double> tri1(a, b, c);
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Delta2::common::triangle(tri1, V, F);
    std::shared_ptr<Delta2::MeshData> M(new Delta2::MeshData(V, F, Delta2::globals::opt, true));

    Eigen::Vector3d a_rot({0, 0, 0});
    Eigen::Vector3d a_pos({0, 2.1, 0});
    Eigen::Matrix4d a_t = Delta2::common::transformationMatrix(a_rot, a_pos);

    Eigen::Vector3d b_rot({0, 0, 0});
    Eigen::Vector3d b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4d b_t = Delta2::common::transformationMatrix(b_rot, b_pos);

    Delta2::Particle PA(M, 1.0, 10.0, 0.25);
    Delta2::Particle PB(M, 1.0, 10.0, 0.25);

    std::vector<Delta2::collision::DeferredCompare> bucket_pairs;
    bucket_pairs.push_back(Delta2::collision::DeferredCompare(M->getSurrogateTree().getNode(0), M->getSurrogateTree().getNode(0)));

    std::vector<Delta2::collision::Contact<double>> hits;
    std::vector<int> pair_used_out;
    Delta2::collision::findContactsBucketHybridDeferred(bucket_pairs, PA, PB, a_t, b_t, hits, pair_used_out);

    REQUIRE( hits.size() == 1 );
    REQUIRE_THAT( hits[0].A, VectorEqual(Eigen::Vector3d({0, 0.7667, 0})).margin(0.01) );
    REQUIRE_THAT( hits[0].B, VectorEqual(Eigen::Vector3d({0, 0.6667, 0})).margin(0.01) );
    REQUIRE( hits[0].eps_a == PA.geo_eps );
    REQUIRE( hits[0].eps_b == PB.geo_eps );
    REQUIRE( hits[0].eps_inner_a == 0 );
    REQUIRE( hits[0].eps_inner_b == 0 );
    REQUIRE( pair_used_out[0] == 1 );
}

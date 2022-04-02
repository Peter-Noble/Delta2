#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "fast_dist.h"
#include "../common/primitive_geo.h"

TEST_CASE("Triangle-triangle fast vert-face", "[grp_collision_detection][fast]") {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    Eigen::Vector3f a_rot({0.5235, 0.3491, 0.0});
    Eigen::Vector3f a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);

    Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};
    
    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_fastInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector3f({-0.342, -0.732, 0.326})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector3f({-0.339, -0.767, 0.129})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector3f({-0.342, -0.732, 0.326}) - Eigen::Vector3f({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
}

TEST_CASE("Triangle-triangle fast edge-edge", "[grp_collision_detection][fast]") {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    Eigen::Vector3f a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Vector3f a_pos({-1.0284, -1.0942, 1.1695});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);

    Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};
    
    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_fastInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector3f({-0.6545, -0.9378, 0.2656})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector3f({-0.5089, -0.8695, 0.1449})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector3f({-0.6545, -0.9378, 0.2656}) - Eigen::Vector3f({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
}

TEMPLATE_TEST_CASE("Triangle-triangle fast intersecting", "[grp_collision_detection][fast]", float, double) {
    Eigen::Vector3f a({-1, 0, 0});
    Eigen::Vector3f b({1, 0, 0});
    Eigen::Vector3f c({0, -2, 0});
    Delta2::common::Triangle<float> tri1(a, b, c);

    Eigen::Vector3f p_out;
    Eigen::Vector3f q_out;

    Eigen::Vector3f a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Vector3f a_pos({0.052572, -0.153374, 1.0443});
    Eigen::Matrix4f a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<float> A = tri1.transformed(a_t);

    Eigen::Vector3f b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector3f b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix4f b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<float> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<float>> Bs = {B};
    
    std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result = Delta2::collision::_fastInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
}

TEST_CASE("Triangle-triangle fast bucket", "[grp_collision_detection][fast]") {
    Delta2::common::Options opt;

    Eigen::MatrixXd V = Eigen::MatrixXd({
        { 1.0,  1.0,  0.0},
        {-1.0,  1.0,  0.0},
        {-1.0, -1.0,  0.0},
        { 1.0, -1.0,  0.0},
        { 0.0,  0.0,  0.1}
    });
    Eigen::MatrixXi F = Eigen::MatrixXi({
        {0, 1, 4},
        {1, 2, 4},
        {2, 3, 4},
        {3, 0, 4},
    });

    std::shared_ptr<Delta2::MeshData> m_surface(new Delta2::MeshData(V, F, opt, true));
    Delta2::Particle p_a(m_surface, 1.0, 1.0, 0.01);
    p_a.future_state.setTranslation({0.0, 0.0, 1.0});
    p_a.is_static = true;

    Delta2::common::plane(1.0, V, F);
    std::shared_ptr<Delta2::MeshData> m_plane(new Delta2::MeshData(V, F, opt, true));
    Delta2::Particle p_b(m_plane, 1.0, 1.0, 0.01);
    p_b.future_state.setTranslation({0.0, 0.0, 1.105});
    
    Delta2::collision::DeferredCompare d(p_a.mesh->getSurrogateTree().getNode(0), p_b.mesh->getSurrogateTree().getNode(0));
    std::vector<Delta2::collision::DeferredCompare> pairs = {d};

    std::vector<Delta2::collision::Contact<double>> hits;

    Delta2::collision::findContactsBucketFast<8, double>(pairs, p_a, p_b, hits);

    REQUIRE( Delta2::collision::filterContacts(hits, 0.0).size() == 1 );
}

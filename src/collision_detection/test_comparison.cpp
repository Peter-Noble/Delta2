#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "comparison_dist.h"
#include "../model/surrogate_components.h"
#include "../common/primitive_geo.h"
#include "../common/viewer.h"

TEMPLATE_TEST_CASE("Triangle-triangle comparison vert-face", "[grp_collision_detection][comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    // Eigen::Vector<TestType, 3> a_rot({0.5235, 0.3491, 0.0}); // Has a problem with rotation order
    Eigen::Quaternion<TestType> a_rot({0.951261, 0.254839, 0.1677497, -0.0449395});
    Eigen::Vector<TestType, 3> a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    // Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternion<TestType> b_rot({0.9952, -0.087056, 0.0038, 0.04347});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_comparisonInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).margin(0.01) );  
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326}) - Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
}

TEMPLATE_TEST_CASE("Triangle-triangle comparison edge-edge", "[grp_collision_detection][comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    // Eigen::Vector<TestType, 3> a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Quaternion<TestType> a_rot({0.744, 0.498, 0.368, 0.252});
    Eigen::Vector<TestType, 3> a_pos({-1.0284, -1.0942, 1.1695});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    // Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternion<TestType> b_rot({0.9952, -0.087056, -0.0038, 0.04347});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_comparisonInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656}) - Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
}

TEMPLATE_TEST_CASE("Triangle-triangle comparison intersecting", "[grp_collision_detection][comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    // Eigen::Vector<TestType, 3> a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Quaternion<TestType> a_rot({0.744, 0.498, 0.368, 0.252});
    Eigen::Vector<TestType, 3> a_pos({0.052572, -0.153374, 1.0443});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    // Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Quaternion<TestType> b_rot({0.9952, -0.087056, -0.0038, 0.04347});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_comparisonInnerDists(A, Bs);

    REQUIRE( result.size() == 1 );
    REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
}

TEST_CASE("Triangle-Triangle comparison bucket", "[grp_collision_detection][comparison]") {
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
    p_a.future_state.setTranslation(Eigen::Vector3d({0.0, 0.0, 1.0}) + m_surface->init_centre_of_mass_offset);
    p_a.is_static = true;

    Delta2::common::plane(1.0, V, F);
    std::shared_ptr<Delta2::MeshData> m_plane(new Delta2::MeshData(V, F, opt, true));
    Delta2::Particle p_b(m_plane, 1.0, 1.0, 0.01);
    p_b.future_state.setTranslation({0.0, 0.0, 1.105});
    
    Delta2::collision::DeferredCompare d(p_a.mesh->getSurrogateTree().getNode(0), p_b.mesh->getSurrogateTree().getNode(0));
    std::vector<Delta2::collision::DeferredCompare> pairs = {d};

    std::vector<Delta2::collision::Contact<double>> hits;
    std::vector<int> pair_used;

    Delta2::collision::findContactsBucketComparison<8, double>(pairs, p_a, p_b, hits, pair_used);

    REQUIRE( pair_used.size() == 1 );
    REQUIRE( pair_used[0] == 8 );
    REQUIRE( Delta2::collision::filterContacts(hits, 0.0).size() == 1 );
}

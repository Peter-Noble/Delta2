#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "penalty_dist.h"

TEMPLATE_TEST_CASE("Triangle-triangle hybrid vert-face", "[grp_collision_detection][hybrid]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    Eigen::Vector<TestType, 3> a_rot({0.5235, 0.3491, 0.0});
    Eigen::Vector<TestType, 3> a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };

    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 8>(A, Bs, max_error, failed);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326}) - Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
}

TEMPLATE_TEST_CASE("Triangle-triangle hybrid edge-edge", "[grp_collision_detection][hybrid]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    Eigen::Vector<TestType, 3> a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Vector<TestType, 3> a_pos({-1.0284, -1.0942, 1.1695});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 8>(A, Bs, max_error, failed);

    REQUIRE( result.size() == 1 );
    REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656})).margin(0.01) );
    REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).margin(0.01) );
    REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656}) - Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
}

TEMPLATE_TEST_CASE("Triangle-triangle hybrid intersecting", "[grp_collision_detection][hybrid]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    Eigen::Vector<TestType, 3> a_rot({1.323341, 0.301492, 0.887539});
    Eigen::Vector<TestType, 3> a_pos({0.052572, -0.153374, 1.0443});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 8>(A, Bs, max_error, failed);

    REQUIRE( result.size() == 1 );
    REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
}

TEMPLATE_TEST_CASE("Triangle-triangle hybrid bucket", "[grp_collision_detection][hybrid]", float, double) {
    REQUIRE( false );
    // TODO test that result is correct within tolerance or correctly identified as not converged
}

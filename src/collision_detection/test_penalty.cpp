#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "penalty_dist.h"
#include "comparison_dist.h"

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-N vert-face", "[grp_collision_detection][penalty]", float, double) {
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

    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 16>(A, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326})).margin(0.01) );
        REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).margin(0.01) );
        REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326}) - Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-N edge-edge", "[grp_collision_detection][penalty]", float, double) {
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
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 16>(A, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656})).margin(0.1) );
        REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).margin(0.1) );
        REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656}) - Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-N intersecting", "[grp_collision_detection][penalty]", float, double) {
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
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 16>(A, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-N range", "[grp_collision_detection][penalty]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;

    auto rx = GENERATE(0.0, 0.17, 0.34, 0.51, 0.68, 0.85, 1.02, 1.19, 1.36);
    auto ry = GENERATE(0.0, 0.17, 0.34, 0.51, 0.68, 0.85, 1.02, 1.19, 1.36);
    auto rz = GENERATE(0.0, 0.17, 0.34, 0.51, 0.68, 0.85, 1.02, 1.19, 1.36);
    Eigen::Vector<TestType, 3> a_rot({rx, ry, rz});
    Eigen::Vector<TestType, 3> a_pos({0.052572, -0.153374, 1.0443});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);
    Delta2::common::Triangle<TestType> A = tri1.transformed(a_t);

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.005 };
    std::vector<bool> failed = { false };
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penaltyInnerDists<TestType, 16>(A, Bs, max_error, failed);
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result_true = Delta2::collision::_comparisonInnerDists(A, Bs);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        if (std::get<0>(result_true[0]) == Approx(0.0).margin(0.0001)) {
            REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.2) );
        } else {
            INFO( "Rotation: " << rx << ", " << ry << ", " << rz );
            REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(std::get<1>(result_true[0])).margin(1.5) );
            REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(std::get<2>(result_true[0])).margin(1.5) );
            REQUIRE( std::get<0>(result[0]) == Approx(std::get<0>(result_true[0])).margin(0.2) );
        }
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-1xN vert-face", "[grp_collision_detection][penalty]", float, double) {
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
    std::vector<Delta2::common::Triangle<TestType>> As = {A};

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };

    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penalty1To1Dists<TestType, 16>(As, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326})).margin(0.01) );
        REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).margin(0.01) );
        REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.342, -0.732, 0.326}) - Eigen::Vector<TestType, 3>({-0.339, -0.767, 0.129})).norm()).margin(0.02) );
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-1xN edge-edge", "[grp_collision_detection][penalty]", float, double) {
    
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
    std::vector<Delta2::common::Triangle<TestType>> As = {A};

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penalty1To1Dists<TestType, 16>(As, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE_THAT( std::get<1>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656})).margin(0.1) );
        REQUIRE_THAT( std::get<2>(result[0]), VectorEqual(Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).margin(0.1) );
        REQUIRE( std::get<0>(result[0]) == Approx((Eigen::Vector<TestType, 3>({-0.6545, -0.9378, 0.2656}) - Eigen::Vector<TestType, 3>({-0.5089, -0.8695, 0.1449})).norm()).margin(0.02) );
    }
}

TEMPLATE_TEST_CASE("Triangle-triangle penalty 1-1xN intersecting", "[grp_collision_detection][penalty]", float, double) {
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
    std::vector<Delta2::common::Triangle<TestType>> As = {A};

    Eigen::Vector<TestType, 3> b_rot({-0.1745, 0.0, 0.0873});
    Eigen::Vector<TestType, 3> b_pos({0.0, 0.0, 0.0});
    Eigen::Matrix<TestType, 4, 4> b_t = Delta2::common::transformationMatrix(b_rot, b_pos);
    Delta2::common::Triangle<TestType> B = tri1.transformed(b_t);
    std::vector<Delta2::common::Triangle<TestType>> Bs = {B};
    
    std::vector<TestType> max_error = { 0.01 };
    std::vector<bool> failed = { false };
    
    std::vector<std::tuple<TestType, Eigen::Vector<TestType, 3>, Eigen::Vector<TestType, 3>>> result = Delta2::collision::_penalty1To1Dists<TestType, 16>(As, Bs, max_error, failed);

    if (!failed[0]) {
        REQUIRE( result.size() == 1 );
        REQUIRE( std::get<0>(result[0]) == Approx(0.0).margin(0.0001) );
    }
}

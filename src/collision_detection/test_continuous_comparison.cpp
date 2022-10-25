#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "continuous_comparison.h"
#include "../common/primitive_geo.h"

TEMPLATE_TEST_CASE("Sphere-sphere first intersect", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a_pos = {0, 0, 0};
    Eigen::Vector<TestType, 3> a_vel = {1, 0, 0};
    TestType a_rad = 1.0;
    Eigen::Vector<TestType, 3> b_pos = {0, -0.5, 0};
    Eigen::Vector<TestType, 3> b_vel = {-1, 0, 0};
    TestType b_rad = 1.0;

    REQUIRE(Delta2::collision::sphereSphereFirstIntersect(a_pos, a_vel, a_rad, b_pos, b_vel, b_rad) == 0.0);

    a_pos = {0, 0, 0};
    a_vel = {1, 0, 0};
    a_rad = 0.5;
    b_pos = {0, -1.001, 0};
    b_vel = {-1, 0, 0};
    b_rad = 0.5;

    REQUIRE(Delta2::collision::sphereSphereFirstIntersect(a_pos, a_vel, a_rad, b_pos, b_vel, b_rad) == -1.0);

    a_pos = {0, 0, 0};
    a_vel = {1, 0, 0};
    a_rad = 0.5;
    b_pos = {0, -1.001, 0};
    b_vel = {1, 0.1, 0};
    b_rad = 0.5;
    REQUIRE(Delta2::collision::sphereSphereFirstIntersect(a_pos, a_vel, a_rad, b_pos, b_vel, b_rad) == Approx(TestType(0.01)).epsilon(0.0001) );
}

TEMPLATE_TEST_CASE("Continuous comparison no intersect far", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;
    TestType toc;
    
    Eigen::Matrix<TestType, 4, 4> a_t_start;
    a_t_start << 0.9961947202682495, -0.0858316496014595, -0.015134435147047043, 0.0,
                 0.08715573698282242, 0.981060266494751, 0.17298738658428192, 0.0,
                 0.0, -0.1736481785774231, 0.9848077297210693, 0.0,
                 0.0, 0.0, 0.0, 1.0;


    Eigen::Matrix<TestType, 4, 4> a_t_end;
    a_t_end << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 2.0,
               0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_start;
    b_t_start << 0.912834644317627, -0.24342688918113708, -0.32783564925193787, -0.8726963400840759,
                 0.40832462906837463, 0.548050582408905, 0.7300079464912415, 1.8028504848480225,
                 0.0019669521134346724, -0.8002399206161499, 0.5996767282485962, 4.524521350860596,
                 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_end;
    b_t_end << 0.9537079334259033, 0.2651883065700531, -0.14183202385902405, 0.0,
               -0.05004240944981575, 0.6049837470054626, 0.7946637272834778, 1.0,
               0.29654160141944885, -0.7507795095443726, 0.5902484655380249, 3.9175784587860107,
               0.0, 0.0, 0.0, 1.0;

    TestType dist = Delta2::collision::triTriCCDLinear(tri1, a_t_start, a_t_end, tri1, b_t_start, b_t_end, p_out, q_out, toc);

    REQUIRE( toc == 1.0 );
    REQUIRE_THAT( p_out, VectorEqual(Eigen::Vector<TestType, 3>({-0.951562, -0.004578, 2.0})).margin(0.1) );
    REQUIRE_THAT( q_out, VectorEqual(Eigen::Vector<TestType, 3>({-0.954302, 1.05403, 3.62159})).margin(0.1) );
    REQUIRE( dist == Approx((Eigen::Vector<TestType, 3>({-0.951562, -0.004578, 2.0}) - Eigen::Vector<TestType, 3>({-0.954302, 1.05403, 3.62159})).norm()).margin(0.05) );
}

TEMPLATE_TEST_CASE("Continuous comparison no intersect close", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;
    TestType toc;
    
    Eigen::Matrix<TestType, 4, 4> a_t_start;
    a_t_start << 0.9961947202682495, -0.0858316496014595, -0.015134435147047043, 0.0,
                 0.08715573698282242, 0.981060266494751, 0.17298738658428192, 0.0,
                 0.0, -0.1736481785774231, 0.9848077297210693, 0.0,
                 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> a_t_end;
    a_t_end << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.43169140815734863,
               0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_start;
    b_t_start << 0.912834644317627, -0.24342688918113708, -0.32783564925193787, 1.958876609802246,
                 0.40832462906837463, 0.548050582408905, 0.7300079464912415, -0.9190752506256104,
                 0.0019669521134346724, -0.8002399206161499, 0.5996767282485962, -1.093824863433838,
                 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_end;
    b_t_end << 0.9537079334259033, 0.2651883065700531, -0.14183202385902405, 1.958876609802246,
               -0.05004240944981575, 0.6049837470054626, 0.7946637272834778, -0.9190752506256104,
               0.29654160141944885, -0.7507795095443726, 0.5902484655380249, 1.0489931106567383,
               0.0, 0.0, 0.0, 1.0;

    TestType dist = Delta2::collision::triTriCCDLinear(tri1, a_t_start, a_t_end, tri1, b_t_start, b_t_end, p_out, q_out, toc);

    REQUIRE( toc == Approx(0.79).margin(0.05) );
    REQUIRE_THAT( p_out, VectorEqual(Eigen::Vector<TestType, 3>({0.64, -0.73, 0.38})).margin(0.05) );
    REQUIRE_THAT( q_out, VectorEqual(Eigen::Vector<TestType, 3>({1.0, -0.95, 0.44})).margin(0.05) );
    REQUIRE( dist == Approx((Eigen::Vector<TestType, 3>({0.64, -0.73, 0.38}) - Eigen::Vector<TestType, 3>({1.0, -0.95, 0.44})).norm()).margin(0.05) );
}

TEMPLATE_TEST_CASE("Continuus comparison intersect", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({-1, 0, 0});
    Eigen::Vector<TestType, 3> b({1, 0, 0});
    Eigen::Vector<TestType, 3> c({0, -2, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> p_out;
    Eigen::Vector<TestType, 3> q_out;
    TestType toc;
    
    Eigen::Matrix<TestType, 4, 4> a_t_start;
    a_t_start << 0.9961947202682495, -0.0858316496014595, -0.015134435147047043, 0.0,
                 0.08715573698282242, 0.981060266494751, 0.17298738658428192, 0.0,
                 0.0, -0.1736481785774231, 0.9848077297210693, 0.0,
                 0.0, 0.0, 0.0, 1.0;


    Eigen::Matrix<TestType, 4, 4> a_t_end;
    a_t_end << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 2.0,
               0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_start;
    b_t_start << 0.9396926164627075, 0.1710100769996643, 0.29619812965393066, 0.0,
                 0.0, 0.8660253882408142, -0.5, 1.0,
                 -0.3420201241970062, 0.46984633803367615, 0.813797652721405, 1.2651395797729492,
                 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<TestType, 4, 4> b_t_end;
    b_t_end << 0.9396926164627075, 0.0, 0.3420201241970062, 0.0,
               0.0, 1.0, 0.0, 1.0,
               -0.3420201241970062, 0.0, 0.9396926164627075, 1.5,
               0.0, 0.0, 0.0, 1.0;

    TestType dist = Delta2::collision::triTriCCDLinear(tri1, a_t_start, a_t_end, tri1, b_t_start, b_t_end, p_out, q_out, toc);

    REQUIRE( toc == Approx(0.24).margin(0.05) );
    REQUIRE_THAT( p_out, VectorEqual(Eigen::Vector<TestType, 3>({-0.265017, -0.843225, 0.591865})).margin(0.1) );
    REQUIRE_THAT( q_out, VectorEqual(Eigen::Vector<TestType, 3>({-0.265017, -0.843225, 0.591865})).margin(0.1) );
    REQUIRE( dist == 0.0 );
}


TEMPLATE_TEST_CASE("Edge edge Newton-Raphson intersect", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({1, -1, 0});
    Eigen::Vector<TestType, 3> b({-1, 1, 0});
    Eigen::Vector<TestType, 3> c({-1, -1, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> P;
    Eigen::Vector<TestType, 3> Q;
    TestType toc;
    
    Eigen::Matrix<TestType, 4, 4> a_T_start = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(66.2817), Delta2::common::degToRad(-45.7037), Delta2::common::degToRad(-111.169)}), Eigen::Vector<TestType, 3>({0, 0, 0}));
    Eigen::Matrix<TestType, 4, 4> a_T_end = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(66.2817), Delta2::common::degToRad(-45.7037), Delta2::common::degToRad(-111.169)}), Eigen::Vector<TestType, 3>({0, 0, 0.7}));
    
    Eigen::Matrix<TestType, 4, 4> b_T_start = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(261.326), Delta2::common::degToRad(52.4147), Delta2::common::degToRad(-101.34)}), Eigen::Vector<TestType, 3>({0, 0, 1.3}));
    Eigen::Matrix<TestType, 4, 4> b_T_end = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(261.326), Delta2::common::degToRad(52.4147), Delta2::common::degToRad(-101.34)}), Eigen::Vector<TestType, 3>({0, 0, 0.5}));
    
    Delta2::common::Triangle<TestType> a_start = tri1.transformed(a_T_start);
    Delta2::common::Triangle<TestType> a_end = tri1.transformed(a_T_end);
    Delta2::common::Triangle<TestType> b_start = tri1.transformed(b_T_start);
    Delta2::common::Triangle<TestType> b_end = tri1.transformed(b_T_end);

    TestType dist = Delta2::collision::edgeEdgeCCD(a_start.AB(), a_end.AB(), b_start.AB(), b_end.AB(), 10000.0, P, Q, toc);

    REQUIRE( toc == Approx(0.87).margin(0.01) );
    REQUIRE_THAT( P, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.6})).margin(0.01) );
    REQUIRE_THAT( Q, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.6})).margin(0.01) );
    REQUIRE( dist == Approx(0.0).margin(0.01) );
}

TEMPLATE_TEST_CASE("Edge edge Newton-Raphson no intersect", "[grp_collision_detection][continuous_comparison]", float, double) {
    Eigen::Vector<TestType, 3> a({1, -1, 0});
    Eigen::Vector<TestType, 3> b({-1, 1, 0});
    Eigen::Vector<TestType, 3> c({-1, -1, 0});
    Delta2::common::Triangle<TestType> tri1(a, b, c);

    Eigen::Vector<TestType, 3> P;
    Eigen::Vector<TestType, 3> Q;
    TestType toc;
    
    Eigen::Matrix<TestType, 4, 4> a_T_start = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(66.2817), Delta2::common::degToRad(-45.7037), Delta2::common::degToRad(-111.169)}), Eigen::Vector<TestType, 3>({0, 0, 0}));
    Eigen::Matrix<TestType, 4, 4> a_T_end = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(66.2817), Delta2::common::degToRad(-45.7037), Delta2::common::degToRad(-111.169)}), Eigen::Vector<TestType, 3>({0, 0, 0.5}));
    
    Eigen::Matrix<TestType, 4, 4> b_T_start = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(261.326), Delta2::common::degToRad(52.4147), Delta2::common::degToRad(-101.34)}), Eigen::Vector<TestType, 3>({0, 0, 1.3}));
    Eigen::Matrix<TestType, 4, 4> b_T_end = Delta2::common::transformationMatrix(Eigen::Vector<TestType, 3>({Delta2::common::degToRad(261.326), Delta2::common::degToRad(52.4147), Delta2::common::degToRad(-101.34)}), Eigen::Vector<TestType, 3>({0, 0, 0.7}));
    
    Delta2::common::Triangle<TestType> a_start = tri1.transformed(a_T_start);
    Delta2::common::Triangle<TestType> a_end = tri1.transformed(a_T_end);
    Delta2::common::Triangle<TestType> b_start = tri1.transformed(b_T_start);
    Delta2::common::Triangle<TestType> b_end = tri1.transformed(b_T_end);

    TestType dist = Delta2::collision::edgeEdgeCCD(a_start.AB(), a_end.AB(), b_start.AB(), b_end.AB(), 100000, P, Q, toc);

    REQUIRE( toc == Approx(1.0).margin(0.02) );
    REQUIRE_THAT( P, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.5})).margin(0.1) );
    REQUIRE_THAT( Q, VectorEqual(Eigen::Vector<TestType, 3>({-0.01, -0.11, 0.68})).margin(0.1) );
    REQUIRE( dist == Approx(0.2).margin(0.05) );
}


TEST_CASE("Continuous bucket", "[grp_collision_detection][continuous_comparison]") {
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
    p_a.current_state.setTranslation({0.0, 0.0, 1.0});
    p_a.future_state.setTranslation({0.0, 0.0, 1.0});
    p_a.is_static = true;

    Delta2::common::plane(1.0, V, F);
    std::shared_ptr<Delta2::MeshData> m_plane(new Delta2::MeshData(V, F, opt, true));
    Delta2::Particle p_b(m_plane, 1.0, 1.0, 0.01);
    p_b.current_state.setTranslation({0.0, 0.0, 2.0});
    p_b.future_state.setTranslation({0.0, 0.0, 1.0});
    
    Delta2::collision::DeferredCompare d(p_a.mesh->getSurrogateTree().getNode(0), p_b.mesh->getSurrogateTree().getNode(0));
    std::vector<Delta2::collision::DeferredCompare> pairs = {d};

    std::vector<Delta2::collision::ContinuousContact<double>> hits;
    std::vector<int> num_hits;

    Delta2::collision::findContactsBucketContinuousComparison<8, double>(pairs, p_a, p_b, hits, num_hits);

    std::vector<Delta2::collision::ContinuousContact<double>> hits_filtered = Delta2::collision::filterContinuousContacts(hits, 0.0);

    REQUIRE( hits_filtered.size() == 1 );
    REQUIRE( hits_filtered[0].toc == Approx(0.9) );
}

TEST_CASE("pointTriCCDLinear error case", "[grp_collision_detection][continuous_comparison]") {
    Delta2::common::Edge<double> e({2, -1, 0.29}, {2, -1, -0.81});
    Delta2::common::Triangle<double> tri({50, 50,0}, {-50, 50, 0}, {50, -50, 0});

    bool isIntersecting;
    double t;
    bool t_inf;
    Eigen::Vector<double, 3> intersect = tri.intersectSegment(e, t, t_inf, isIntersecting);

    REQUIRE(isIntersecting);
}

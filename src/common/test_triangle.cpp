#include "triangle.h"
#include "../testing/unit_test.h"
#include "../testing/approx_triangle.h"

TEMPLATE_TEST_CASE( "Triangle init and retrieve", "[grp_common][triangle]", float, double ) {
    Eigen::Vector<TestType, 3> a = {0.0, 1.1, 2.2};
    Eigen::Vector<TestType, 3> b = {3.0, 4.1, 5.2};
    Eigen::Vector<TestType, 3> c = {6.0, 7.1, 8.2};
    Delta2::common::Edge<TestType> ab(a, b);
    Delta2::common::Edge<TestType> bc(b, c);
    Delta2::common::Edge<TestType> ca(c, a);
    Delta2::common::Triangle<TestType> T(a, b, c);
    REQUIRE( T.A == a );
    REQUIRE( T.B == b );
    REQUIRE( T.C == c );

    REQUIRE( T[0] == a );
    REQUIRE( T[1] == b );
    REQUIRE( T[2] == c );
    REQUIRE_THROWS( T[-1] );
    REQUIRE_THROWS( T[3] );

    REQUIRE( T.AB().A == a );
    REQUIRE( T.AB().B == b );
    REQUIRE( T.BC().A == b );
    REQUIRE( T.BC().B == c );
    REQUIRE( T.CA().A == c );
    REQUIRE( T.CA().B == a );
}

TEMPLATE_TEST_CASE( "Triangle property methods", "[grp_common][triangle]", float, double ) {
    Delta2::common::Triangle<TestType> pt({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    REQUIRE( pt.volume(1.0) == Approx(4.0/3.0 * 3.14159) );
    REQUIRE( pt.area() == 0.0 );
    REQUIRE( pt.diameter() == 0.0 );

    Eigen::Vector<TestType, 3> a({0.0, 0.0, 0.0});
    Eigen::Vector<TestType, 3> b({1.0, 0.0, 0.0});
    Eigen::Vector<TestType, 3> c({0.0, 1.0, 0.0});
    Delta2::common::Triangle<TestType> tri(a, b, c);

    REQUIRE( tri.area() == Approx(0.5) );
    REQUIRE_THAT( tri.normal(), VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 1.0})) );
    REQUIRE( tri.diameter() == Approx(1.414213562) );

    REQUIRE_THAT( tri.template cast<float>(), TriangleEqual(Delta2::common::Triangle<float>(a.template cast<float>(), b.template cast<float>(), c.template cast<float>())) );
    REQUIRE_THAT( tri.template cast<double>(), TriangleEqual(Delta2::common::Triangle<double>(a.template cast<double>(), b.template cast<double>(), c.template cast<double>())) );
}


TEMPLATE_TEST_CASE( "Triangle interaction methods", "[grp_common][triangle]", float, double ) {
    Eigen::Vector<TestType, 3> a({0.0, 0.0, 0.0});
    Eigen::Vector<TestType, 3> b({1.0, 0.0, 0.0});
    Eigen::Vector<TestType, 3> c({0.0, 1.0, 0.0});
    Delta2::common::Triangle<TestType> tri(a, b, c);

    auto [x0, y0, z0] = tri.calcBarycentric(a);
    auto [x1, y1, z1] = tri.calcBarycentric(b);
    auto [x2, y2, z2] = tri.calcBarycentric(c);
    REQUIRE( x0 == TestType(1.0) );
    REQUIRE( y0 == TestType(0.0) );
    REQUIRE( z0 == TestType(0.0) );
    REQUIRE( x1 == TestType(0.0) );
    REQUIRE( y1 == TestType(1.0) );
    REQUIRE( z1 == TestType(0.0) );
    REQUIRE( x2 == TestType(0.0) );
    REQUIRE( y2 == TestType(0.0) );
    REQUIRE( z2 == TestType(1.0) );

    Eigen::Vector<TestType, 3> d({-1.0, 1.0, 0.0});
    Eigen::Vector<TestType, 3> e({-1.0, -1.0, 0.0});

    Eigen::Vector<TestType, 3> start({-1.0, 1.0, 1.0});
    Eigen::Vector<TestType, 3> end({1.0, -1.0, -1.0});
    
    Delta2::common::Triangle<TestType> tri2(b, d, e);
    Delta2::common::Edge<TestType> E(start, end);
    TestType t;
    bool isIntersecting;
    bool t_inf;
    Eigen::Vector<TestType, 3> intersect_pt = tri2.intersectSegment(E, t, t_inf, isIntersecting);

    REQUIRE( isIntersecting );
    REQUIRE( t == Approx(t) );
    REQUIRE_THAT( intersect_pt, VectorEqual(Eigen::Vector<TestType, 3>({0, 0, 0})) );

    Eigen::Vector<TestType, 3> p({2.0, 2.0, 1.0});
    Eigen::Vector<TestType, 3> p_on_plane({2.0, 2.0, 0.0});

    Eigen::Vector<TestType, 3> project_pt = tri2.projectToPlane(p, isIntersecting);

    REQUIRE( !isIntersecting );
    REQUIRE_THAT( project_pt, VectorEqual(p_on_plane) );

    Eigen::Vector<TestType, 3> in_v0({2.0, 0, 0});

    Delta2::common::TriangleZone zone = tri2.projectPointToZone(in_v0);

    REQUIRE( zone == Delta2::common::TriangleZone::V0 );

    TestType dist = tri2.distToPoint(in_v0);

    REQUIRE( dist == Approx(1.0) );
}

#include "edge.h"
#include "../testing/unit_test.h"
#include "../testing/approx_edge.h"

TEMPLATE_TEST_CASE( "Edge init and retrieve", "[grp_common][edge]", float, double ) {
    Eigen::Vector<TestType, 3> a = {0.0, 1.1, 2.2};
    Eigen::Vector<TestType, 3> b = {3.0, 4.1, 5.2};
    Delta2::common::Edge<TestType> ab(a, b);
    REQUIRE( ab.A == a );
    REQUIRE( ab.B == b );

    Eigen::Vector<float, 3> af = {0.0, 1.1, 2.2};
    Eigen::Vector<float, 3> bf = {3.0, 4.1, 5.2};

    Eigen::Vector<double, 3> ad = {0.0, 1.1, 2.2};
    Eigen::Vector<double, 3> bd = {3.0, 4.1, 5.2};

    REQUIRE_THAT( ab.template cast<float>(), EdgeEqual(Delta2::common::Edge<float>(af, bf)) );
    REQUIRE_THAT( ab.template cast<double>(), EdgeEqual(Delta2::common::Edge<double>(ad, bd)) );
}

TEMPLATE_TEST_CASE( "Edge interaction methods", "[grp_common][edge]", float, double ) {
    Eigen::Vector<TestType, 3> a = {0.0, 1.0, 0.0};
    Eigen::Vector<TestType, 3> b = {0.0, -1.0, 0.0};
    Delta2::common::Edge<TestType> ab(a, b);

    Eigen::Vector<TestType, 3> c = {1.0, 0.0, 1.0};
    Eigen::Vector<TestType, 3> d = {1.0, 0.0, -1.0};
    Delta2::common::Edge<TestType> cd(c, d);
    
    Eigen::Vector<TestType, 3> e = {0.0, 0.0, 1.0};
    Eigen::Vector<TestType, 3> f = {0.0, 0.0, 0.0};
    Delta2::common::Edge<TestType> ef(e, f);

    REQUIRE( ab.dist(cd) == Approx(1.0) );
    REQUIRE( ab.dist(ef) == Approx(0.0) );
    REQUIRE( ef.dist(cd) == Approx(1.0) );

    REQUIRE( ab.dist(cd) == Approx(cd.dist(ab)) );
    REQUIRE( ab.dist(ef) == Approx(ef.dist(ab)) );
    REQUIRE( ef.dist(cd) == Approx(cd.dist(ef)) );

    Eigen::Vector<TestType, 3> p_out, q_out;
    Eigen::Vector<TestType, 3> p_out_sym, q_out_sym;
    ab.closest(cd, p_out, q_out);
    cd.closest(ab, p_out_sym, q_out_sym);
    REQUIRE_THAT( p_out, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.0})) );
    REQUIRE_THAT( q_out, VectorEqual(Eigen::Vector<TestType, 3>({1.0, 0.0, 0.0})) );
    REQUIRE_THAT( p_out, VectorEqual(q_out_sym) );
    REQUIRE_THAT( q_out, VectorEqual(p_out_sym) );

    ab.closest(ef, p_out, q_out);
    ef.closest(ab, p_out_sym, q_out_sym);
    REQUIRE_THAT( p_out, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.0})) );
    REQUIRE_THAT( q_out, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 0.0, 0.0})) );
    REQUIRE_THAT( p_out, VectorEqual(q_out_sym) );
    REQUIRE_THAT( q_out, VectorEqual(p_out_sym) );
    
    ef.closest(cd, p_out, q_out);
    cd.closest(ef, p_out_sym, q_out_sym);
    REQUIRE( p_out.x() == Approx(0.0) );
    REQUIRE( p_out.y() == Approx(0.0) );
    REQUIRE( p_out.z() >= Approx(0.0) );
    REQUIRE( p_out.z() <= Approx(1.0) );

    REQUIRE( q_out.x() == Approx(1.0) );
    REQUIRE( q_out.y() == Approx(0.0) );
    REQUIRE( q_out.z() >= Approx(0.0) );
    REQUIRE( q_out.z() <= Approx(1.0) );

    REQUIRE_THAT( p_out, VectorEqual(q_out_sym) );
    REQUIRE_THAT( q_out, VectorEqual(p_out_sym) );
}

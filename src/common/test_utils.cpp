#include "basic_utils.h"
#include "../testing/unit_test.h"
#include "../testing/approx_matrix.h"
#include "../testing/approx_vector.h"
#include "../common/triangle.h"
#include "utils.h"

TEMPLATE_TEST_CASE("lerp", "[grp_common][basic_utils]", float, double) {
    Eigen::Vector<TestType, 3> A = {0, 0, 0};
    Eigen::Vector<TestType, 3> B = {1, 1, 1};
    
    Eigen::Vector<TestType, 3> r0 = Delta2::common::lerp(A, B, TestType(0.0));
    Eigen::Vector<TestType, 3> r1 = Delta2::common::lerp(A, B, TestType(1.0));
    Eigen::Vector<TestType, 3> r2 = Delta2::common::lerp(A, B, TestType(0.5));

    REQUIRE( r0 == A );
    REQUIRE( r1 == B );
    REQUIRE( r2 == Eigen::Vector<TestType, 3>({TestType(0.5), TestType(0.5), TestType(0.5)}) );
}

TEMPLATE_TEST_CASE("clamp01", "[grp_common][basic_utils]", float, double) {
    REQUIRE( Delta2::common::clamp01(TestType(-1)) == TestType(0.0));
    REQUIRE( Delta2::common::clamp01(TestType(-0.5)) == TestType(0.0));
    REQUIRE( Delta2::common::clamp01(TestType(0)) == TestType(0.0));
    REQUIRE( Delta2::common::clamp01(TestType(0.5)) == TestType(0.5));
    REQUIRE( Delta2::common::clamp01(TestType(1)) == TestType(1.0));
    REQUIRE( Delta2::common::clamp01(TestType(1.5)) == TestType(1.0));
    REQUIRE( Delta2::common::clamp01(TestType(2)) == TestType(1.0));
}

TEMPLATE_TEST_CASE("transformationMatrix rotation", "[grp_common][basic_utils]", float, double) {
    Eigen::Vector<TestType, 3> a_rot({0.0, 0.0, 3.14159/2.0});
    Eigen::Vector<TestType, 3> a_pos({0.0, 0.0, 0});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);

    Eigen::Vector<TestType, 3> pt({1.0, 0.0, 0.0});

    Eigen::Matrix<TestType, 4, 4> a_t_true;
    a_t_true << 0.0, -1.0, 0.0, 0.0,
                1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0;
    REQUIRE_THAT( a_t, MatrixEqual(a_t_true) );

    Eigen::Vector<TestType, 3> result = Delta2::common::transform(pt, a_t);

    REQUIRE_THAT( result, VectorEqual(Eigen::Vector<TestType, 3>({0.0, 1.0, 0.0})) );
}

TEMPLATE_TEST_CASE("transformationMatrix rotation 2", "[grp_common][basic_utils]", float, double) {
    // Eigen::Vector<TestType, 3> a_rot({0.5235, 0.3491, 0.0}); // Has a problem with rotation order
    Eigen::Quaternion<TestType> a_rot({0.951261, 0.254839, 0.1677497, -0.0449395});
    Eigen::Vector<TestType, 3> a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);

    Eigen::Vector<TestType, 3> pt({1.0, 0.0, 0.0});

    Eigen::Matrix<TestType, 4, 4> a_t_true;
    a_t_true << 0.9396809339523315, 0.17099685966968536, 0.2962428033351898, 0.0,
                0.0, 0.866074800491333, -0.4999144673347473, 1.0,
                -0.3420522212982178, 0.46976011991500854, 0.8138339519500732, 1.2651000022888184,
                0.0, 0.0, 0.0, 1.0;
    REQUIRE_THAT( a_t, MatrixEqual(a_t_true) );

    Eigen::Vector<TestType, 3> result = Delta2::common::transform(pt, a_t);

    REQUIRE_THAT( result, VectorEqual(Eigen::Vector<TestType, 3>({0.9397, 1.0, 0.923})) );
}

TEMPLATE_TEST_CASE("transformationMatrix", "[grp_common][basic_utils]", float, double) {
    // Eigen::Vector<TestType, 3> a_rot({0.5235, 0.3491, 0.0}); // Has a problem with rotation order
    Eigen::Quaternion<TestType> a_rot({0.951261, 0.254839, 0.1677497, -0.0449395});
    Eigen::Vector<TestType, 3> a_pos({0.0, 1.0, 1.2651});
    Eigen::Matrix<TestType, 4, 4> a_t = Delta2::common::transformationMatrix(a_rot, a_pos);

    Eigen::Matrix<TestType, 4, 4> a_t_true;
    a_t_true << 0.9396926164627075, 0.1710100769996643, 0.29619812965393066, 0.0,
                0.0, 0.8660253882408142, -0.5, 1.0,
                -0.3420201241970062, 0.46984633803367615, 0.813797652721405, 1.2651395797729492,
                0.0, 0.0, 0.0, 1.0;
    REQUIRE_THAT( a_t, MatrixEqual(a_t_true) );
}

TEST_CASE("toTriangles", "[grp_common][utils]") {
    Eigen::MatrixXd V {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 0.0}
    };
    Eigen::MatrixXi F {
        {0, 1, 2},
        {3, 1, 2}
    };
    std::vector<Delta2::common::Triangle<double>> tris = Delta2::common::toTriangles<double>(V, F, F.rows());
    REQUIRE( tris.size() == 2 );
    REQUIRE( tris[0].A == Eigen::Vector3d(V(0,0), V(0,1), V(0,2)) );
    REQUIRE( tris[0].B == Eigen::Vector3d(V(1,0), V(1,1), V(1,2)) );
    REQUIRE( tris[0].C == Eigen::Vector3d(V(2,0), V(2,1), V(2,2)) );

    REQUIRE( tris[1].A == Eigen::Vector3d(V(3,0), V(3,1), V(3,2)) );
    REQUIRE( tris[1].B == Eigen::Vector3d(V(1,0), V(1,1), V(1,2)) );
    REQUIRE( tris[1].C == Eigen::Vector3d(V(2,0), V(2,1), V(2,2)) );

    REQUIRE( tris[0].area() + tris[1].area() == Approx(1.0) );
}


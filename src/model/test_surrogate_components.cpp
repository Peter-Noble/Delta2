
// TEMPLATE_TEST_CASE("BucketConnected basic - size 2", "[grp_model][bucket]", float, double) {
//     Eigen::MatrixXd V {
//         {0.0, 0.0, 0.0},
//         {1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         {1.0, 1.0, 0.0}
//     };
//     Eigen::MatrixXi F {
//         {0, 1, 2},
//         {3, 1, 2}
//     };

//     Delta2::model::BucketConnected bucket(V.template cast<TestType>(), F);
//     std::vector<Delta2::common::Triangle<TestType>> tris = bucket.getTriangles();

//     REQUIRE( tris.size() == 2 );
//     REQUIRE( bucket.getEdges().size() == 5 );
// }

// TEMPLATE_TEST_CASE("BucketSoup basic - size 2", "[grp_model][bucket]", float, double) {
//     Eigen::MatrixXd V {
//         {0.0, 0.0, 0.0},
//         {1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         {1.0, 1.0, 0.0}
//     };
//     Eigen::MatrixXi F {
//         {0, 1, 2},
//         {3, 1, 2}
//     };

//     Delta2::model::BucketSoup<TestType, 2> bucket(V.template cast<TestType>(), F);
//     std::vector<Delta2::common::Triangle<TestType>> tris = bucket.getTriangles();

//     REQUIRE( tris.size() == 2 );
//     REQUIRE( bucket.getEdges().size() == 6 );
// }

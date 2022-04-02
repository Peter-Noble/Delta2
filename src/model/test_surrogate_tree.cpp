#include "../testing/unit_test.h"
#include "surrogate_tree.h"
#include "../common/primitive_geo.h"

// TEMPLATE_TEST_CASE("Surrogate tree - single bucket", "[grp_model][surrogate_tree]", float, double) {
//     Delta2::common::Options opt;

//     Eigen::MatrixXd Vd;
//     Eigen::MatrixXi F;

//     Delta2::common::plane(1.0, Vd, F);
//     Eigen::Matrix<TestType, -1, -1> V = Vd.cast<TestType>();

//     Delta2::model::SurrogateTree st(V, F, opt);

//     REQUIRE( st.getNumNodes() == 1 );
//     std::shared_ptr<Delta2::model::Bucket> b = st.getBucket(0);
//     REQUIRE( b->getTriangles().size() == 2 );
//     REQUIRE( b->getEps().size() == 2 );
//     REQUIRE( b->getEps()[0] == 0.0 );
//     REQUIRE( b->getEps()[1] == 0.0 );
// }

// TEMPLATE_TEST_CASE("Surrogate tree - sphere", "[grp_model][surrogate_tree]", float, double) {
//     Delta2::common::Options opt;

//     Eigen::MatrixXd Vd;
//     Eigen::MatrixXi F;

//     Delta2::common::sphere(1.0, 5, Vd, F);
//     Eigen::Matrix<TestType, -1, -1> V = Vd.cast<TestType>();

//     Delta2::model::SurrogateTreeTemplate<TestType, 8, 8> st(V, F, opt);

//     REQUIRE( st.getNumNodes() == 9 );
//     const Delta2::model::Node& n = st.getNode(0);
//     REQUIRE( n.bucket_id == 0 );
//     REQUIRE( n.depth == 0 );
//     REQUIRE( n.is_inner );
//     REQUIRE( n.num_children == 8 );
//     REQUIRE( n.parent_id == -1 );

//     std::shared_ptr<Delta2::model::Bucket> b = st.getBucket(0);
//     REQUIRE( b->getTriangles().size() == 8 );
//     for (const Delta2::common::Triangle<TestType>& t : b->getTriangles()) {
//         REQUIRE( t.area() > 0.0 );
//     }
//     for (int i = 1; i <= 8; i++) {
//         const Delta2::model::Node& child_n = st.getNode(i);
//         std::shared_ptr<Delta2::model::Bucket> child_b = st.getBucket(i);
//         REQUIRE( child_n.parent_id == 0 );
//         REQUIRE( child_n.depth == 1 );
//         std::shared_ptr<Delta2::model::BucketSoup<TestType, 8>> bucket = std::dynamic_pointer_cast<Delta2::model::BucketSoup<TestType, 8>>(b);
//         const Delta2::model::BucketSoup<TestType, 8>& bucket_unpacked = *bucket.get();
//         if (child_b->getTriangles().size() > 1) {
//             REQUIRE( b->getEps()[i-1] > 0.0 );
//         }
//     }
// }

// TEMPLATE_TEST_CASE("Surrogate tree - sphere large", "[grp_model][surrogate_tree]", float, double) {
//     Delta2::common::Options opt;

//     Eigen::MatrixXd Vd;
//     Eigen::MatrixXi F;

//     Delta2::common::sphere(1.0, 16, Vd, F);
//     Eigen::Matrix<TestType, -1, -1> V = Vd.cast<TestType>();

//     Delta2::model::SurrogateTreeTemplate<TestType, 8, 8> st(V, F, opt);

//     REQUIRE( st.getNumNodes() > (16-2)*16/8 );
// }

// TEMPLATE_TEST_CASE("Surrogate tree - qslim sphere large", "[grp_model][surrogate_tree]", float, double) {
//     REQUIRE( false );
// }

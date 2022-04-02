#include "../testing/unit_test.h"
#include "../testing/approx_vector.h"
#include "convex_hull.h"
#include <stdlib.h>


TEST_CASE("Convex hull 1", "[grp_common][convex_hull]") {
    std::vector<Eigen::Vector3d> pts;
    pts.push_back(Eigen::Vector3d({0.707733, 0.593858, 0.382683}));
    pts.push_back(Eigen::Vector3d({0.173648, 0.984808, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.16043, 0.909844, 0.382683}));
    pts.push_back(Eigen::Vector3d({0.707733, 0.593858, 0.382683}));
    pts.push_back(Eigen::Vector3d({0.766044, 0.642788, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.173648, 0.984808, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.766044, 0.642788, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.16043, 0.909844, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.173648, 0.984808, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.766044, 0.642788, 6.12323e-17}));
    pts.push_back(Eigen::Vector3d({0.707733, 0.593858, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.16043, 0.909844, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.707733, 0.593858, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.122788, 0.696364, -0.707107}));
    pts.push_back(Eigen::Vector3d({0.16043, 0.909844, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.707733, 0.593858, -0.382683}));
    pts.push_back(Eigen::Vector3d({0.541675, 0.454519, -0.707107}));
    pts.push_back(Eigen::Vector3d({0.122788, 0.696364, -0.707107}));
    
    std::vector<Delta2::common::Triangle<double>> hull = Delta2::common::convex_hull(pts);

    REQUIRE(hull.size() > 0);
}

TEST_CASE("Convex hull 2", "[grp_common][convex_hull]") {
    std::vector<Eigen::Vector3d> pts;

    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < 2; y++) {
            for (int z = 0; z < 2; z++) {
                pts.push_back(Eigen::Vector3d({x, y, z}));
            }
        }
    }

    for (int i = 0; i < 100; i++) {
        pts.push_back(Eigen::Vector3d({float(rand())/float((RAND_MAX)), float(rand())/float((RAND_MAX)), float(rand())/float((RAND_MAX))}));
    }

    std::vector<Delta2::common::Triangle<double>> hull = Delta2::common::convex_hull(pts);

    Delta2::common::bbox bbox = hull[0].bbox();

    for (Delta2::common::Triangle<double>& t : hull) {
        bbox = bbox.expand(t.bbox());
    }

    VectorApproxEqual<double> lower(Eigen::Vector3d({0.0, 0.0, 0.0}));
    VectorApproxEqual<double> upper(Eigen::Vector3d({1.0, 1.0, 1.0}));

    REQUIRE(lower.match(bbox.lower));
    REQUIRE(upper.match(bbox.upper));
}
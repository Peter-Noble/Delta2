#include "../testing/unit_test.h"
#include "particle.h"
#include "../common/primitive_geo.h"
#include "../testing/approx_matrix.h"

TEST_CASE("Inertia matrix", "[grp_model][mesh_data]") {
    Delta2::common::Options opt;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::Matrix3d I, I_true;

    Delta2::common::cube(1, V, F);
    Delta2::MeshData M(V, F, opt);

    I = M.getUnitInertiaTensor();
    double cube_const = 1.0/6.0 * 4;
    I_true = M.getVolume() * Eigen::Matrix3d({
        {cube_const,          0,          0},
        {         0, cube_const,          0},
        {         0,          0, cube_const}
    });

    CHECK_THAT( I, MatrixEqual(I_true) );

    Delta2::common::sphere(1.0, 10, V, F);
    M = Delta2::MeshData(V, F, opt);

    I = M.getUnitInertiaTensor();
    double sphere_const = 2.0/5.0;
    I_true = M.getVolume() * Eigen::Matrix3d({
        {sphere_const,            0,            0},
        {           0, sphere_const,            0},
        {           0,            0, sphere_const}
    });

    CHECK_THAT( I, MatrixEqual(I_true).margin(0.15) );

    Delta2::common::sphere(1.0, 8, V, F);
    M = Delta2::MeshData(V, F, opt, true);

    I = M.getUnitInertiaTensor();
    double hollow_sphere_const = 2.0/3.0;
    I_true = M.getVolume() * Eigen::Matrix3d({
        {hollow_sphere_const,                   0,                   0},
        {                  0, hollow_sphere_const,                   0},
        {                  0,                   0, hollow_sphere_const}
    });

    CHECK_THAT( I, MatrixEqual(I_true).margin(0.15) );
}
